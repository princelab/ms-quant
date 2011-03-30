#!/usr/bin/env ruby

require 'ms/ident/peptide_hit/qvalue'
require 'ms/ident/protein_hit'
require 'ms/ident/peptide/db'
require 'ms/quant/spectral_counts'

require 'trollop'

def putsv(*args)
  if $VERBOSE
    puts(*args) ; $stdout.flush
  end
end

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} peptide_centric_db.yml, file1.psq ...
}
  opt :names, "array of names for the table (otherwise filenames)", :type => String
  opt :fdr_percent, "%FDR as cutoff", :default => 1.0
  opt :write_subset, "(development) write subset db", :default => false
end

opt = opts.parse(ARGV)

if ARGV.size < 2
  opts.educate && exit
end

peptide_centric_db_file = ARGV.shift

opt[:names] ||= ARGV.map do |file| 
  base = file.chomp(File.extname(file))
  base=base.chomp(File.extname(base)) if File.extname(base) == '.phq'
  base
end

class Ms::Ident::PeptideHit
  attr_accessor :experiment_name
end
fdr_cutoff = opt[:fdr_percent] / 100

start=Time.now

$VERBOSE = true

ar_of_peptide_hit_ars = Ms::Ident::Peptide::Db::IO.open(peptide_centric_db_file) do |peptide_to_proteins|
  putsv "#{Time.now-start} seconds to read #{peptide_centric_db_file}"
  ARGV.zip(opt[:names]).map do |file,exp| 
    peptide_hits = Ms::Ident::PeptideHit::Qvalue.from_file(file)
    putsv "#{file}: #{peptide_hits.size} hits"
    peptide_hits.select! do |hit|
      if hit.qvalue <= fdr_cutoff
        # update each peptide with its protein hits
        prot_ids = peptide_to_proteins[hit.aaseq]
        if prot_ids
          hit.experiment_name = exp
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
  ar_of_peptide_hit_ars.each do |pephits|
    pephits.each do |pephit|
      aaseqs_to_prots[pephit.aaseq] = pephit.proteins
    end
  end
  outfile = "peptidecentric_subset.yml"
  puts "writing #{outfile} with #{aaseqs_to_prots.size} aaseq->protids"
  File.open(outfile,'w') do |out| 
    aaseqs_to_prots.each do |k,v| 
      out.puts(%Q{#{k}: #{v.map(&:id).join("\t") }}) 
    end
  end
end

$VERBOSE = true
if $VERBOSE
  opt[:names].zip(ar_of_peptide_hit_ars) do |name, pep_ar|
    puts "#{name}: #{pep_ar.size}"
  end
end

all_peptide_hits = ar_of_peptide_hit_ars.flatten(1)


# because peptide_hit#proteins yields id strings (which hash properly),
# each protein group is an array of 
protein_groups = Ms::Ident::ProteinGroup.peptide_hits_to_protein_groups(all_peptide_hits)

pephit_to_protein_groups = Hash.new {|h,k| h[k] = [] }
protein_groups.each do |protein_group|
  protein_group.peptide_hits.each {|hit| pephit_to_protein_groups[hit] << protein_group }
end

# partition them all out by filename

ar_of_count_data = opt[:names].map do |name|
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

# protein_groups
# [ ar_of_counts_for_exp1, ar_of_counts_for_exp2, ar_of_counts_for_exp3 ]

protein_groups.zip(*ar_of_count_data) do |row|
  pg = row.shift
  puts (row.map(&:to_a).flatten + pg.to_a).join("\t")
end


