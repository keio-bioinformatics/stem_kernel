#!/usr/bin/env ruby
require 'dishuffle'
require 'normal'

$avg_rate=0.25
$div_rate=0.05
$lim_rate=0.5

def rand_len(l)
  avg=l*$avg_rate
  div=l*$div_rate
  x=p_norm(rand-0.5)*div+avg
  x=0 if x<0.0
  x=l*$lim_rate if x>l*$lim_rate
  x.to_i
end  

gets('>')
while s=gets('>')
  l=s.split("\n")
  desc=l.shift
  s=l.join("")
  s=s[0..-2] if s[-1..-1]==">"
  up_r=dishuffle(s)
  down_r=dishuffle(s)
  up_len=rand_len(s.size)
  down_len=rand_len(s.size)
  seq=up_r[up_r.size/2,up_len]+s+down_r[down_r.size/2,down_len]
  puts ">#{desc} (orig #{s.size}, upstream #{up_len}, downstream #{down_len}, total #{up_len+down_len+s.size})",seq
end
