#!/usr/bin/env ruby

require 'dishuffle'

con_rank=2

seq=Hash.new("")
name=[]

# load an alignment
while l=gets
  next if l=~/multiple sequence alignment/
  next if l=~/^\s+$/
  n,s=l.chomp.split
  next if s=~/[^ACGTUacgtu\-]/
  name.push(n) unless seq.include?(n)
  seq[n]+=s
end

len=0
ma=[]
name.each do |n|
  s=seq[n]
  len=s.size if len==0
  raise if len!=s.size
  s.split('').each_with_index{|v,i|
    ma[i]=[] if ma[i].nil?
    ma[i].push(v)
  }
end

th=(name.size.to_f*0.5).to_i
r = ma.map do |v|
  hist=Hash.new(0)
  v.each{|vv| hist[vv]+=1}
  vv=hist.keys.sort{|a,b| hist[b]<=>hist[a] }
  if hist[vv[0]]>th
    vv[0]
  else
    vv[0,con_rank].join('')
  end
end

q,idx=dishuffle(r)

puts "CLUSTAL W (1.83) multiple sequence alignment"
puts

while idx.size>0
  x=nil
  if idx.size>50
    x=idx[0,50]
    idx=idx[50..-1]
  else
    x=idx
    idx=[]
  end

  puts
  name.each_with_index do |n,j|
    puts n.ljust(25)+x.map{|i| ma[i][j]}.join('')
  end
  puts
end
