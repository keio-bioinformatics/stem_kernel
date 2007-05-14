#!/usr/bin/env ruby

require 'optparse'

n_seqs=0

opt = OptionParser.new
opt.on('-n seqs', 'the number of seqs in each output file'){|v| n_seqs=v.to_i}
opt.parse!(ARGV)

seq=[]
gets(">")
while fa=gets(">")
  seq.push(fa.gsub(/>$/m,''))
end

n_seqs = seq.size if seq.size<n_seqs
(0..n_seqs-1).each do |i|
  print ">"+seq[i]
end
