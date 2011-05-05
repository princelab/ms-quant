#!/usr/bin/env ruby

require 'andand'
require 'set'
require 'ruport'

require 'ms/ident/peptide_hit/qvalue'
require 'ms/ident/protein'
require 'ms/ident/peptide/db'
require 'ms/quant/spectral_counts'
require 'ms/quant/protein_group_comparison'
require 'ms/quant/qspec/protein_group_comparison'
require 'ms/quant/qspec'
require 'ms/quant/cmdline'


require 'yaml'
require 'tempfile'

require 'trollop'

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

class Ruport::Data::Table
  # returns self
  def add_column_with_data(colname, array_of_data, opts={})
    self.add_column(colname, opts)
    self.data.zip(array_of_data) do |row, newval|
      row[colname] = newval
    end
    self
  end

  # acceptable opts:
  #
  #     :header => an array of lines (each which will be commented out)
  def to_tsv(file, opt={})
    delimiter = "\t" 
    File.open(file,'w') do |out|
      opt[:header].each {|line| out.puts "# #{line}" } if opt[:header]
      out.puts self.column_names.join(delimiter)
      self.data.each do |row|
        out.puts row.to_a.join(delimiter)
      end
      opt[:footer].each {|line| out.puts "# #{line}" } if opt[:footer]
    end
  end

end

def write_subset(sample_to_pephits, outfile="peptidecentric_subset.yml")
  aaseqs_to_prots = {} 
  sample_to_pephits.map(&:last).flatten(1).each do |pephit|
    aaseqs_to_prots[pephit.aaseq] = pephit.proteins.map(&:id)
  end
  File.open(outfile,'w') do |out| 
    aaseqs_to_prots.each do |k,v| 
      out.puts(%Q{#{k}: #{v.join("\t") }}) 
    end
  end
end


outfile = "spectral_counts.tsv"
pephits_outfile = "spectral_counts.pephits.tsv"
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
  opt :peptides, "also write peptide hits (to: #{pephits_outfile})", :default => false
  opt :verbose, "speak up", :default => false
  opt :count_type, "type of spectral counts (<spectral|aaseqcharge|aaseq>)", :default => 'spectral'
  opt :qspec_normalize, "normalize spectral counts per run", :default => false
  opt :qspec_keep_files, "keep a copy of the files submitted and returned from Qspec", :default => false
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

(samplename_to_filename, condition_to_samplenames, samplename_to_condition) = Ms::Quant::Cmdline.args_to_hashes(ARGV)

raise ArgumentError, "must have 2 conditions for qspec!" if opt[:qspec] && condition_to_samplenames.size != 2

samplenames = samplename_to_filename.keys

class Ms::Ident::PeptideHit
  attr_accessor :experiment_name
  attr_accessor :protein_groups
end
class Ms::Ident::Protein
  attr_accessor :length
end


fdr_cutoff = opt[:fdr_percent] / 100

if opt[:qspec] || opt[:descriptions]
  putsv "reading lengths and descriptions from #{opt[:fasta]}"
  (id_to_length, id_to_desc) = Ms::Fasta.protein_lengths_and_descriptions(opt[:fasta])
end

samplename_to_peptidehits = samplename_to_filename.map do |sample, file|
 [sample, Ms::Ident::PeptideHit::Qvalue.from_file(file).select {|hit| hit.qvalue <= fdr_cutoff }]
end

# update each peptide hit with protein hits and sample name:
all_protein_hits = Hash.new {|h,id| h[id] = Ms::Ident::Protein.new(id) }
Ms::Ident::Peptide::Db::IO.open(peptide_centric_db_file) do |peptide_to_proteins|
  samplename_to_peptidehits.map do |sample, peptide_hits|
    peptide_hits.each do |hit|
      # update each peptide with its protein hits
      protein_hits = peptide_to_proteins[hit.aaseq].map do |id| 
        protein = all_protein_hits[id]
        protein.length = id_to_length[id] if id_to_length
        protein.description = id_to_desc[id] if id_to_desc
        protein
      end
      hit.experiment_name = sample
      # if there are protein hits, the peptide hit is selected 
      hit.proteins = protein_hits
    end
  end
end

write_subset(samplename_to_peptidehits) if opt[:write_subset]

samplename_to_peptidehits.each {|samplename, hits| putsv "#{samplename}: #{hits.size}" } if $VERBOSE

all_peptide_hits = samplename_to_peptidehits.map(&:last).flatten(1)

# this constricts everything down to a minimal set of protein groups that
# explain the entire set of peptide hits.   
update_pephits = true  # ensures that each pephit is linked to the array of protein groups it is associated with
protein_groups = Ms::Ident::ProteinGroup.peptide_hits_to_protein_groups(all_peptide_hits, update_pephits)

hits_table_hash = {}  # create the table using key => column hash
samplenames.each do |name|
  hits_table_hash[name] = protein_groups.map do |prot_group|
    prot_group.peptide_hits.select {|hit| hit.experiment_name == name }
  end
end

# The columns are filled with groups of peptide hits, one group of hits per
# protein group (protein group order is implicit).  The rows are sample names.
#
#  (implied) sample1   sample2   sample3   ...
#  (group1)  [hit,hit] [hit...]  [hit...]  ...
#  (group2)  [hit,hit] [hit...]  [hit...]  ...
#   ...         ...       ...       ...       ...
hits_table = Ruport::Data::Table.new(:data => hits_table_hash.values.transpose, :column_names => hits_table_hash.keys)

# spectral counts of type opt[:count_type]
counts_data = hits_table.data.map do |row|
  row.map do |pephits|
    Ms::Quant::SpectralCounts.counts(pephits) {|pephit| 1.0 / pephit.protein_groups.size }.send(opt[:count_type])
  end
end

# each cell holds a SpectralCounts object, which hash 3 types of count data
counts_table = Ruport::Data::Table.new(:data => counts_data, :column_names => samplenames)

# return a list of ProteinGroupComparisons
if opt[:qspec]

  # prepare data for qspec
  condition_to_count_array = counts_table.column_names.map do |name| 
    [samplename_to_condition[name], counts_table.column(name)] 
  end
  # average length of the proteins in the group
  name_length_pairs = protein_groups.map do |pg|
    [pg.map(&:id).join(":"), pg.map(&:length).reduce(:+)./(pg.size).round]
  end

  qspec_results = Ms::Quant::Qspec.new(name_length_pairs, condition_to_count_array).run(opt[:qspec_normalize], :keep => opt[:qspec_keep_files])
  
  cols_to_add = [:bayes_factor, :fold_change, :fdr]
  counts_table.add_columns cols_to_add
  counts_table.data.zip(qspec_results) do |row, qspec_result|
    cols_to_add.each {|cat| row[cat] = qspec_result[cat] }
  end
end

counts_table.add_columns( [:name, :ids, :description] )
counts_table.data.zip(protein_groups) do |row, pg|
  best_id = pg.sort_by {|prot| [prot.id, prot.length] }.first
  row.name = best_id.description.andand.match(/ GN=([^\s]+) ?/).andand[1] || best_id.id
  row.ids = pg.map(&:id).join(',')
  row.description = best_id.description
end


if opt[:peptides]
  hits_table.each do |record|
    record.each_with_index do |hits,i|
      new_cell = hits.group_by do |hit| 
        [hit.aaseq, hit.charge]
      end.map do |key, hits|
        [key.reverse.join("_"), hits.map(&:id).join(',')].join(":")
      end.join('; ')
      record[i] = new_cell
    end
  end
  hits_table.add_column_with_data(:name, counts_table.column(:name), :position=>0)
  hits_table.to_tsv(pephits_outfile, :footer => ["parallel to #{outfile}"])
end

intro = ["samples: #{samplename_to_filename}", "options: #{opt}"]
counts_table.to_tsv(outfile, :footer => intro)
