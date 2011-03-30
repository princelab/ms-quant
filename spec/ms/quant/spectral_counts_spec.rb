require 'spec_helper'

require 'ms/quant/spectral_counts'



PeptideHit = Struct.new(:aaseq, :charge, :proteins) do
  def initialize(*args)
    super(*args)
    self.proteins ||= []
  end
  def inspect # easier to read output
    "<PeptideHit aaseq=#{self.aaseq} charge=#{self.charge} proteins(ids)=#{self.proteins.map(&:id).join(',')}>"
  end
  def hash ; self.object_id end
end
ProteinHit = Struct.new(:id, :peptide_hits) do
  def inspect # easier to read output
    "<Prt #{self.id}>"
  end
  def hash ; self.object_id end
end

describe 'groups of peptide hits' do
  before do
    @pep_hits = [ 
      ['AABBCCDD', 2], #bg,mg1,mg2 0.33
      ['BBCC', 2],     #bg,mg1,mg2 0.33
      ['DDEEFFGG', 2], #bg,mg1,mg2,sbm 0.25
      ['DDEEFFGG', 2], #bg,mg1,mg2,sbm 0.25
      ['DDEEFFGG', 3], #bg,mg1,mg2,sbm 0.25
      ['HIYA', 2],     #bg,lg 0.5
    ].map {|ar| PeptideHit.new(*ar) }
    @prot_hits = {
                              # spectral_counts, aaseq+charge counts, aaseq counts
      'big_guy' => @pep_hits, # 6, 5, 4; 
      'little_guy' => [@pep_hits.last], # 1, 1, 1, 0.5, 0.5, 0.5
      'medium_guy1' => @pep_hits[0,5], # 5, 4, 3 
      'medium_guy2' => @pep_hits[0,5], # 5, 4, 3 
      'subsumed_by_medium' => @pep_hits[2,3], # 3, 2, 1
    }.map {|data| ProteinHit.new(*data) }
    # doubly linked for this
    @prot_hits.each do |prot| 
      prot.peptide_hits.each {|pephit| pephit.proteins << prot }
    end
    # DEPENDS ON AN ORDERED HASH (RUBY 1.9!!!!)
    @expected_counts = [ [6,5,4], [1,1,1], [5,4,3], [5,4,3], [3,2,1] ]
    @expected_counts_split = [ [1.9167,1.6667,1.4167], [0.5,0.5,0.5], [1.41667, 1.1667, 0.91667], [1.41667, 1.16667, 0.91667], [0.75, 0.5, 0.25] ]
  end

  it 'finds spectral counts (without sharing)' do
    groups_of_pephits = @prot_hits.map(&:peptide_hits)
    counts = Ms::Quant::SpectralCounts.counts(groups_of_pephits)
    @expected_counts.zip(counts) do |exp, act|
      act.to_a.is exp
    end
  end

  it 'finds spectral counts (splitting counts between shared)' do
    groups_of_pephits = @prot_hits.map(&:peptide_hits)
    counts = Ms::Quant::SpectralCounts.counts(groups_of_pephits) {|pephit| pephit.proteins.size }
    @expected_counts_split.zip(counts) do |exp, act|
      exp.zip(act) {|e,a| a.should.be.close e, 0.0001 }
    end
  end

end
