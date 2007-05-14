#!/usr/bin/env ruby

require 'roc'

fold=[]

while l=gets
  if l=~/^== (\d+) ([+-]?\d+) ([+-]?[\d.]+)/
    fl=$1.to_i
    fold[fl]=[] if fold[fl].nil?
    fold[fl].push([$2.to_i,$3.to_f])
  elsif l=~/^Cross/
    print l
  end
end

sum=0.0
sum2=0.0
num=0
fold.each do |f|
  next if f.nil?
  s,=roc(f)
  sum+=s*f.size
  sum2+=s*s*f.size
  num+=f.size
end
avg=sum/num
var=sum2/num-avg*avg
var=0.0 if var<0.0
puts "ROC score = #{avg}, #{Math.sqrt(var)}"
