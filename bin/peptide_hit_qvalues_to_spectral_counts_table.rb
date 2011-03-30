#!/usr/bin/env ruby

require 'ms/ident/peptide_hit/qvalue'
require 'ms/ident/protein_hit'
require 'ms/ident/peptide/db'
require 'ms/quant/spectral_counts'

require 'trollop'

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} peptide_centric_db.yml, file1.psq ...
}
  opt :names, "array of names for the table (otherwise filenames)", :type => String
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

start=Time.now
ar_of_peptide_hit_ars = Ms::Ident::Peptide::Db::IO.open(peptide_centric_db_file) do |peptide_to_proteins|
  puts "#{Time.now-start} seconds to read #{peptide_centric_db_file}" ; $stdout.flush

  hit_to_exp = {}
  ARGV.zip(opt[:names]).map do |file,exp| 
    peptide_hits = Ms::Ident::PeptideHit::Qvalue.from_file(file)
    peptide_hits.select! do |hit|
      # update each peptide with its protein hits
      hit.proteins = peptide_to_proteins[hit].map {|id| Ms::Ident::ProteinHit.new(id) }
    end
    peptide_hits.each {|h| hit_to_exp[hit] = exp }
    peptide_hits
  end
end

all_peptide_hits = ar_of_peptide_hit_ars.flatten(1)
protein_groups = Ms::Ident::ProteinGroup.peptide_hits_to_protein_groups(all_peptide_hits)

ar_of_count_data = opt[:names].map do |name|
  pep_hit_to_prot_groups = Hash.new {|h,k| h[k] = [] }
  groups_of_pephits = protein_groups.map do |prot_group|
    pep_hits = peptide_hits.select {|hit| hit_to_exp[hit] == name }
    pep_hits.each do |pep_hit|
      pep_hit_to_prot_groups[pep_hit] << prot_group
    end # returns the group of pep_hits
  end
  counts = Ms::Quant::SpectralCounts.counts(groups_of_pephits, pep_hit_to_prot_groups)
end

# protein_groups
# [ ar_of_counts_for_exp1, ar_of_counts_for_exp2, ar_of_counts_for_exp3 ]

protein_groups.zip(*ar_of_count_data) do |row|
  p row
end


