#!/usr/bin/env ruby

require 'ms/ident/peptide_hit/qvalue'
require 'ms/ident/protein_hit'
require 'ms/ident/peptide/db'
require 'ms/quant/spectral_counts'
require 'ms/quant/qspec'

require 'yaml'
require 'tempfile'

require 'trollop'

# inverse from Tilo Sloboda (now in facets)

class Hash
  def inverse
    i = Hash.new
    self.each_pair do |k,v|
      if (Array === v) ; v.each{ |x| i[x] = ( i.has_key?(x) ? [k,i[x]].flatten : k ) }
      else ; i[v] = ( i.has_key?(v) ? [k,i[v]].flatten : k ) end
    end ; i  
  end
end


def putsv(*args)
  if $VERBOSE
    puts(*args) ; $stdout.flush
  end
end

def basename(file)
  base = file.chomp(File.extname(file))
  base=base.chomp(File.extname(base)) if File.extname(base) == '.phq'
  base
end


outfile = "spectral_counts.tsv"
delimiter = "\t"

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} <fasta>.peptide_centric_db.yml group1=f1.psq,f2.psq group2=f3.psq,f4.psq
or (each file a group):    #{File.basename(__FILE__)} <fasta>.peptide_centric_db.yml file1.psq file2.psq ...

writes to #{outfile}
group names can be arbitrarily defined
psq is really .psq.tsv file
}
  opt :fdr_percent, "%FDR as cutoff", :default => 1.0
  opt :qspec, "return qspec results (executes qspec or qspecgp). Requires :fasta.  Only 2 groups currently allowed", :default => false
  opt :descriptions, "include descriptions of proteins, requires :fasta", :default => false
  opt :fasta, "the fasta file.  Required for :qspec and :descriptions", :type => String
  opt :outfile, "the to which file data are written", :default => outfile
  opt :verbose, "speak up", :default => false
  opt :count_type, "type of spectral counts (<spectral|aaseqcharge|aaseq>)", :default => 'spectral'
  opt :qspec_normalize, "normalize spectral counts per run", :default => false
  opt :write_subset, "(dev use only) write subset db", :default => false
end

opt = opts.parse(ARGV)
opt[:count_type] = opt[:count_type].to_sym

$VERBOSE = opt.delete(:verbose)

if ARGV.size < 2
  opts.educate && exit
end

if (opt[:qspec] || opt[:descriptions]) && !opt[:fasta]
  puts "You must provide a fasta file with --fasta to use qspec or descriptions!!"
  opts.educate && exit
end

peptide_centric_db_file = ARGV.shift
raise ArgumentError, "need .yml file for peptide centric db" unless File.extname(peptide_centric_db_file) == '.yml'
putsv "using: #{peptide_centric_db_file} as peptide centric db"

# groupname => files
condition_to_samplenames = {}
samplename_to_filename = {}
ARGV.each do |arg|
  (condition, files) = 
    if arg.include?('=')
      (condition, filestring) = arg.split('=')
      [condition, filestring.split(',')]
    else
      [basename(arg), [arg]]
    end
  reptag = ARGV.size
  sample_to_file_pairs = files.each_with_index.map {|file,i| ["#{condition}-rep#{i+1}", file] }
  sample_to_file_pairs.each {|name,file| samplename_to_filename[name] = file }
  condition_to_samplenames[condition] = sample_to_file_pairs.map(&:first)
end


if $VERBOSE
  puts "** condition: sample_names"
  puts condition_to_samplenames.to_yaml
  puts "** samplename: filename"
  puts samplename_to_filename.to_yaml
end

raise ArgumentError, "must have 2 conditions for qspec!" if opt[:qspec] && condition_to_samplenames.size != 2

samplenames = samplename_to_filename.keys

class Ms::Ident::PeptideHit
  attr_accessor :experiment_name
end
fdr_cutoff = opt[:fdr_percent] / 100

start=Time.now

ar_of_pephit_ars = Ms::Ident::Peptide::Db::IO.open(peptide_centric_db_file) do |peptide_to_proteins|
  putsv "#{Time.now-start} seconds to read #{peptide_centric_db_file}"
  samplename_to_filename.map do |sample, file|
    peptide_hits = Ms::Ident::PeptideHit::Qvalue.from_file(file)
    putsv "#{file}: #{peptide_hits.size} hits"
    peptide_hits.select! do |hit|
      if hit.qvalue <= fdr_cutoff
        # update each peptide with its protein hits
        prot_ids = peptide_to_proteins[hit.aaseq]
        if prot_ids
          hit.experiment_name = sample
          hit.proteins = prot_ids
        else ; false end
      else
        false
      end
    end
    peptide_hits
  end
end

if opt[:write_subset]
  aaseqs_to_prots = {} 
  ar_of_pephit_ars.flatten(1).each do |pephit|
    aaseqs_to_prots[pephit.aaseq] = pephit.proteins
  end
  outfile = "peptidecentric_subset.yml"
  puts "writing #{outfile} with #{aaseqs_to_prots.size} aaseq->protids"
  File.open(outfile,'w') do |out| 
    aaseqs_to_prots.each do |k,v| 
      out.puts(%Q{#{k}: #{v.join("\t") }}) 
    end
  end
end

if $VERBOSE
  samplenames.zip(ar_of_pephit_ars) do |samplename, pep_ar|
    putsv "#{samplename}: #{pep_ar.size}"
  end
end

all_peptide_hits = ar_of_pephit_ars.flatten(1)

# because peptide_hit#proteins yields id strings (which hash properly),
# each protein group is an array of 
protein_groups = Ms::Ident::ProteinGroup.peptide_hits_to_protein_groups(all_peptide_hits)

pephit_to_protein_groups = Hash.new {|h,k| h[k] = [] }
protein_groups.each do |protein_group|
  protein_group.peptide_hits.each {|hit| pephit_to_protein_groups[hit] << protein_group }
end

# partition them all out by filename

counts_parallel_to_names_with_counts_per_group = samplenames.map do |name|
  pep_hit_to_prot_groups = Hash.new {|h,k| h[k] = [] }
  groups_of_pephits = protein_groups.map do |prot_group|
    pep_hits = prot_group.peptide_hits.select {|hit| hit.experiment_name == name }
    pep_hits.each do |pep_hit|
      pep_hit_to_prot_groups[pep_hit] << prot_group
    end # returns the group of pep_hits
  end
  counts = Ms::Quant::SpectralCounts.counts(groups_of_pephits) # do |pephit|
  #  pephit_to_protein_groups[pephit].size
  #end
end

if opt[:qspec] || opt[:descriptions]
  putsv "reading lengths and descriptions from #{opt[:fasta]}"
  (id_to_length, id_to_desc) = Ms::Fasta.protein_lengths_and_descriptions(opt[:fasta])
end

samplename_to_condition = condition_to_samplenames.inverse

### OUTPUT TABLE
header_cats = samplenames.map.to_a

ar_of_rows = counts_parallel_to_names_with_counts_per_group.map do |counts_per_group|
  counts_per_group.map(&opt[:count_type])
end.transpose

if opt[:qspec]
  all_conditions = samplenames.map {|sn| samplename_to_condition[sn] }
  condition_to_count_array = all_conditions.zip(counts_parallel_to_names_with_counts_per_group).map do |condition, counts_par_groups|
    [condition, counts_par_groups.map(&opt[:count_type])]
  end

  name_length_pairs = protein_groups.map do |pg|
    # prefer swissprot (sp) proteins over tremble (tr) and shorter protein
    # lengths over longer lengths
    best_guess_protein_id = pg.sort_by {|prot_id| [prot_id, -id_to_length[prot_id]] }.first
    length = id_to_length[best_guess_protein_id]
    [pg.join(":"), length]
  end

  putsv "qspec to normalize counts: #{opt[:qspec_normalize]}"
  qspec_results = Ms::Quant::Qspec.new(name_length_pairs, condition_to_count_array).run(opt[:qspec_normalize])

  to_add = [:fdr, :bayes_factor, :fold_change]
  header_cats.push(*to_add)
  qspec_results.zip(ar_of_rows) do |zipped|
    (result, row) = zipped
    row.push(*to_add.map {|v| result.send(v) })
  end
end

header_cats.push( *%w(BestID AllIDs) )
header_cats.push( 'Description' ) if opt[:descriptions]

protein_groups.zip(ar_of_rows) do |zipped|
  (pg, row) = zipped
  # swiss-prot and then the shortest
  best_protid = pg.sort_by {|prot_id| [prot_id, -id_to_length[prot_id]] }.first
  (gene_id, desc) = 
    if opt[:descriptions] 
      desc = id_to_desc[best_protid] 
      gene_id = (md=desc.match(/ GN=(\w+) ?/)) ? md[1] : best_protid
      [gene_id, desc]
    else
      [best_protid, nil]
    end
  row << gene_id << pg.join(',')
  row.push(desc) if desc
end

### SORT???

File.open(opt[:outfile],'w') do |out|
  out.puts header_cats.join(delimiter) 
  ar_of_rows.each {|row| out.puts row.join(delimiter) }
  putsv "wrote: #{opt[:outfile]}"
end

