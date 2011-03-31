module Ms ; end
module Ms::Quant ; end

class Ms::Quant::Qspec

  NBURNIN = 50  # need to check
  NITER = 2000 # check
  INIT_HEADER = %w(protid protLen)
  DELIMITER = "\t"

  # takes an ordered list of conditions ['cond1', 'cond1', 'cond2', 'cond2'] and
  # returns an array of ints [0,0,0,1,1,1...]
  def self.conditions_to_ints(conditions)
    i = 0
    current_condition = conditions.first
    conditions.map do |cond| 
      if current_condition == cond ; i
      else
        i += 1
        current_condition = cond
        i
      end
    end
  end

  # returns an array of Results structs which is each row of the returned file
  # works with V2 of QSpec
  def self.results_array(resultsfile)
    rows = IO.readlines(resultsfile).map {|line| line.chomp.split("\t") }
    headers = rows.shift
    rows.map do |row| 
      data = [row[0]]
      data.push( *row[1,2].map(&:to_i) )
      data.push( *row[3,4].map(&:to_f) )
      data.push( row[7] )
      Results.new(*data)
    end
  end

  # returns the right executable based on the array of conditions
  def self.executable(conditions)
    biggest_size = conditions.group_by {|v| v }.values.map(&:size).max
    (biggest_size >= 3) ? 'qspecgp' : 'qspec'
  end

  # protname_length_pairs is an array of doublets: [protname, length]
  # condition_to_count_array is an array doublets: [condition, array_of_counts]
  def initialize(protname_length_pairs, condition_to_count_array)
    @protname_length_pairs = protname_length_pairs
    @condition_to_count_array = condition_to_count_array
  end

  def conditions
    @condition_to_count_array.map(&:first)
  end

  # writes a qspec formatted file to filename
  def write(filename)
    ints = Ms::Quant::Qspec.conditions_to_ints(conditions)
    header_cats = INIT_HEADER + ints
    rows = @protname_length_pairs.map {|pair| pair.map.to_a }
    @condition_to_count_array.each do |cond,counts|
      rows.zip(counts) {|row,cnt| row << cnt }
    end
    File.open(filename,'w') do |out|
      out.puts header_cats.join(DELIMITER)
      rows.each {|row| out.puts row.join(DELIMITER) }
    end
  end

  def run(normalize=true, opts={})
    puts "normalize: #{normalize}" if $VERBOSE
    tfile = Tempfile.new("qspec")
    write(tfile.path)
    qspec_exe = self.class.executable(conditions)
    cmd = [qspec_exe, tfile.path, NBURNIN, NITER, (normalize ? 1 : 0)].join(' ')
    puts "running #{cmd}" if $VERBOSE
    reply = `#{cmd}`
    puts reply if $VERBOSE
    results = self.class.results_array(tfile.path + Results::EXT)
    tfile.unlink
    results
  end

  # for version 2 of QSpec
  Results = Struct.new(:protid, :set0, :set1, :bayes_factor, :fold_change, :rb_stat, :fdr, :flag)
  class Results
    EXT = '_qspec'
  end
end
