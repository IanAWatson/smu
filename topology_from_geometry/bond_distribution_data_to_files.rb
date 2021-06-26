# frozen_string_literal: true

# Convert the raw output from the bond length distribution
# to individual files.
# input looks like:

# 6.0.6.13282.3
# 6.0.6.15806.1
# 8.0.6.19522.1
# 6.0.6.15781.1
# 6.0.6.19942.1
# 6.1.1.5865.31
# 6.1.1.5957.16
# 6.1.1.5795.12
# 6.0.1.11350.10
# 6.0.1.11327.8

# Group by the first three tokkens (atype.btype.atype) and then
# for each bonded pair, generate the bond length distribution.
# The 4th token is the distance * 10.000

require 'docopt'

doc = <<~DOCOPT
    Convert output of get_bond_length_distribution to individual files
  #{'  '}
    Usage:
      #{__FILE__} --input=<input> --stem=<stem>
  #{'  '}
    Options:
      --input=<input>     Specify input file
      --stem=<stem>       File name stem for data files generated
  #{'  '}
DOCOPT

begin
  require 'pp'
  pp Docopt.docopt(doc)
rescue Docopt::Exit => e
  puts e.message
end

options = Docopt.docopt(doc)

# Return the index of the first non-nil and last non-nil entries in `data`
def first_and_last(data)
  first_non_nil = -1
  last_non_nil = -1
  data.each_with_index do |value, ndx|
    next unless value
    next unless ndx > 10_000

    first_non_nil = ndx unless first_non_nil >= 0
    last_non_nil = ndx
  end

  [first_non_nil, last_non_nil]
end

# Write a single bond length distribution to `fname`.
# Args:
#  fname: the output file
#  data: a vector of counts. We only write from the first to last
#        non zero values.
#        The indices are converted by distances by dividing by 10,000
def write_bond_length_distribution(fname, data) # rubocop:disable Metrics/MethodLength
  first_ndx, last_ndx = first_and_last(data)
  if (last_ndx - first_ndx) < 2
    $stderr << "Ignoring no data #{fname}\n"
    return
  end
  $stderr << "Opening #{fname} range #{first_ndx} to #{last_ndx}\n"
  File.open(fname, 'w') do |output|
    (first_ndx..last_ndx).each do |ndx|
      output << ndx.to_f / 10_000.0
      output << if data[ndx]
                  ",#{data[ndx]}\n"
                else
                  ",0\n"
                end
    end
  end
end

def main(input_fname, stem) # rubocop:disable Metrics/MethodLength, Metrics/AbcSize
  bonded_pair = {}
  File.open(input_fname, 'r').each do |line|
    f = line.chomp.split('.')
    bpair = f[0..2].join('.')
    unless bonded_pair.key?(bpair)
      $stderr << "Created pair for #{bpair}\n"
      bonded_pair[bpair] = []
    end
    bonded_pair[bpair][f[3].to_i] = f[4]
  end

  $stderr << "Read #{bonded_pair.size} bonded pair data\n"

  bonded_pair.each do |key, value|
    fname = "#{stem}.#{key}"
    write_bond_length_distribution(fname, value)
  end
end

main(options['--input'], options['--stem'])
