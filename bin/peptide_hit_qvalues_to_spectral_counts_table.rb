#!/usr/bin/env ruby

require 'ms/ident/peptide_hit/qvalue'
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

opt[:names] ||= ARGV.map {|file| file.chomp(File.extname(file)) }

peptide_to_proteins = Ms::Ident::Peptide::Db::IO.open(peptide_centric_db_file)

hit_to_exp = {}
ar_of_peptide_hit_ars = ARGV.zip(opt[:names]).map do |file,exp| 
  peptide_hits = Ms::Ident::PeptideHit::Qvalue.from_file(file)
  peptide_hits.select! do |hit|
    # update each peptide with its protein hits
    hit.proteins = peptide_to_proteins[hit]
  end
  peptide_hits.each {|h| hit_to_exp[hit] = exp }
  peptide_hits
end

all_peptide_hits = ar_of_peptide_hit_ars.flatten(1)
protein_groups = Ms::Ident::ProteinGroup.peptide_hits_to_protein_groups(all_peptide_hits)

# prot_group -> { label -> peptide_hit_array }
hit_to_prot_group = Hash.new {|h,k| h[k] = [] }
protein_groups.each do |prot_group|
  prot_group.peptide_hits.group_by {|hit| hit_to_exp[hit] }
end



Ms::Quant::SpectralCounts.counts


