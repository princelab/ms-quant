
module Ms
  module Quant
  end
end

module Ms::Quant::ProteinGroupComparison

  # a protein group object
  attr_accessor :protein_group

  # an array of experiment names
  attr_accessor :experiments

  # parallel array to experiments with the measured values
  attr_accessor :values

  def initialize(protein_group, experiments, values)
    (@protein_group, @experiment, @values) = protein_group, experiment, values
  end
end

class Ms::Quant::ProteinGroupComparison::SpectralCounts
  include Ms::Quant::ProteinGroupComparison
end

class Ms::Quant::ProteinGroupComparison::UniqAAzCounts
  include Ms::Quant::ProteinGroupComparison
end
