require 'set'
require 'ms/ident/protein_group'

module Ms
  module Quant
    module SpectralCounts
      Counts = Struct.new(:spectral, :aaseq, :aaseqcharge)

      SORT_PROTEINS = lambda do |group|
        [group.spectral_counts, group.aaseqcharge_counts, group.aaseq_counts]
      end

      # returns a parallel array of Count objects.  If hit_to_object hash
      # given, then counts are split between groups sharing the hit. 
      # peptide_hits must respond to :charge and :aaseq
      def self.counts(groups_of_pephits, hit_to_object=nil)
        if split_counts
          peptide_hit_to_objects = Hash.new {|h,k| h[k] = [] }
          aaseq_to_objects = Hash.new {|h,k| h[k] = Set.new }
          aaseqcharge_to_objects = Hash.new {|h,k| h[k] = Set.new }
          peptide_hits.each do |peptide_hit|
            object = hit_to_object[peptide_hit]
            peptide_hit_to_objects[peptide_hit] << object
            aaseq_to_objects[peptide_hit.aaseq] << object
            aaseqcharge_to_objects[peptide_hit.aaseq, peptide_hit.charge] << object
          end
        end
        peptide_hit_groups.map do |peptide_hits|
          uniq_aaseq = Set.new 
          uniq_aaseq_charge = Set.new 
          peptide_hits.each do |peptide_hit|
            pephits.each do |pephit|
              uniq_aaseq_charge << [pephit.aaseq, pephit.charge]
              uniq_aaseq << pephit.aaseq
            end
            counts_data = 
              if split_counts
                [
                  peptide_hits.inject(0.0) {|sum,pephit| sum += (1.0 / peptide_hit_to_objects[pephit].size) },
                  uniq_aaseq.inject(0.0) {|sum,aaseq| sum += (1.0 / aaseq_to_objects[aaseq].size) },
                  uniq_aaseq_charge.inject(0.0) {|sum,aaseqcharge| sum += (1.0 / aaseqcharge_to_objects[aaseqcharge].size) }
                ]
              else
                [ peptide_hits.size, uniq_aaseq.size, uniq_aaseq_charge.size]
              end
          end
          Counts.new(*counts_data)
        end
      end

      # takes an ordered hash of experiment names pointing to arrays of
      # [Ms::Ident::ProteinGroup, Count].  Expects that the protein groups are
      # shared objects between the different sets of experiments groups.
      # returns a names array and a hash with protein_group keys and values
      # which are an array of count objects per label.
      def self.flatten_counts_for_experiments(experiment_to_protein_group_count_pairs, &sort_by)
        group_to_counts = Hash.new {|h,k| h[k] = [] }
        labels = []
        experiment_to_protein_group_count_pairs.each do |name, pairs|
          pairs.each do |protein_group, count_obj|
            group_to_counts[protein_group] << count_obj
          end
          labels << name
        end
        [labels, groups_to_counts]
      end

    end
  end
end
