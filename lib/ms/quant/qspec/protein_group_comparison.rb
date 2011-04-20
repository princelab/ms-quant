
module Ms
  module Quant
    module ProteinGroupComparison
    end
  end
end

class Ms::Quant::ProteinGroupComparison::Qspec
  include Ms::Quant::ProteinGroupComparison

  attr_accessor :bayes_factor, :fold_change, :rb_stat, :fdr, :flag

  # the values are the counts array
  def initialize(protein_group, experiments, values, bayes_factor, fold_change, rb_stat, fdr, flag)
    super(protein_group, experiments, values)
    (@bayes_factor, @fold_change, @rb_stat, @fdr, @flag) = bayes_factor, fold_change, rb_stat, fdr, flag
  end
end

