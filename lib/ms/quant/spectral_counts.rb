require 'set'
require 'ms/ident/protein_group'

module Ms
  module Quant
    module SpectralCounts
      Counts = Struct.new(:spectral, :aaseqcharge, :aaseq)

      # returns a parallel array of Count objects.  If split_hits then counts
      # are split between groups sharing the hit.  peptide_hits must respond
      # to :charge and :aaseq.  If split_hits, then each peptide_hit must
      # respond to :linked_to yielding an object with a :size reflective of
      # the number of shared peptide_hits.
      def self.counts(groups_of_peptide_hits, &share_the_pephit)
        groups_of_peptide_hits.map do |peptide_hits|
          uniq_aaseq = {}
          uniq_aaseq_charge = {}
          linked_sizes = peptide_hits.map do |hit|
            linked_to_size = share_the_pephit ? share_the_pephit.call(hit) : 1
            # these guys will end up clobbering themselves, but the
            # linked_to_size should be consistent if the key is the same
            uniq_aaseq_charge[[hit.aaseq, hit.charge]] = linked_to_size
            uniq_aaseq[hit.aaseq] = linked_to_size
            linked_to_size
          end
          counts_data = [linked_sizes, uniq_aaseq_charge.values, uniq_aaseq.values].map do |array|
            share_the_pephit ?  array.inject(0.0) {|sum,size| sum+=(1.0/size) } : array.size
          end
          Counts.new(*counts_data)
        end
      end

    end
  end
end
