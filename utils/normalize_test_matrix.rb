#!/usr/bin/env ruby
require 'zlib'

def train_diagonal(fh)
  diag=[]
  fh.each_line do |l|
    x=l.chomp.split
    lb=x.shift
    a,n=x[0].split(':')
    n=n.to_i
    m,v=x[n].split(':')
    diag[n]=Math.sqrt(v.to_f)
  end
  diag
end

def test_diagonal(fh)
  ret=[]
  while line=fh.gets
    ret.push(Math.sqrt(line.chomp.to_f))
  end
  ret
end

matfile=ARGV.shift
tr_diag=[]
if matfile=~/gz$/
  Zlib::GzipReader.open(matfile) {|fh| tr_diag=train_diagonal(fh) }
else
  open(matfile) {|fh| tr_diag=train_diagonal(fh) }
end

ts_norm=ARGV.shift
ts_diag=[]
open(ts_norm) {|fh| ts_diag=test_diagonal(fh) }

while l=gets
  raise if ts_diag.empty?
  l=l.chomp.split
  xx=ts_diag.shift
  label=l.shift
  head=l.shift
  print "#{label} #{head} "
  l.each do |x|
    i,v=x.split(":")
    #puts [v, i, tr_diag[i.to_i], xx].join(" ")
    v=v.to_f/tr_diag[i.to_i]/xx
    print "#{i}:#{v} "
  end
  puts
end
