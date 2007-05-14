#!/usr/bin/env ruby

lim=ARGV.shift.to_i

mat=[]
lab=[]

while l=gets
  x=l.chomp.split
  la=x.shift
  n=x.shift.split(/:/)
  #next if n[1].to_i>lim
  print "#{la} 0:#{n[1]} "
  for n in x
    n=n.split(/:/)
    next if n[0].to_i>lim
    print "#{n[0]}:#{n[1]} "
  end
  puts
end
