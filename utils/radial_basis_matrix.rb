#!/usr/bin/env ruby
require 'zlib'

def diagonal(fh)
  diag=[]
  fh.each_line do |l|
    x=l.chomp.split
    lb=x.shift
    a,n=x[0].split(':')
    n=n.to_i
    m,v=x[n].split(':')
    diag[n]=v.to_f
  end
  diag
end

def normalize_matrix(gamma, diag, fh, out=$>)
  fh.each_line do |l|
    x=l.chomp.split
    lb=x.shift
    a,n=x.shift.split(':')
    n=n.to_i
    out.print "#{lb} 0:#{n} "
    x.each do |xx|
      m,v=xx.split(':')
      m=m.to_i
      v=v.to_f
      vv=Math.exp(-gamma*(diag[n]+diag[m]-2*v))
      out.print "#{m}:#{vv} "
    end
    out.puts
  end
end

gamma=ARGV.shift.to_f
matfile=ARGV.shift
diag=[]
if matfile=~/gz$/
  Zlib::GzipReader.open(matfile) {|fh| diag=diagonal(fh) }
  Zlib::GzipReader.open(matfile) {|fh| normalize_matrix(gamma, diag, fh) }
else
  open(matfile) {|fh| diag=diagonal(fh) }
  open(matfile) {|fh| normalize_matrix(gamma, diag, fh) }
end
