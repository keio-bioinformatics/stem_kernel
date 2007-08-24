#!/usr/bin/env ruby

require 'dishuffle'

gets('>')
while s=gets('>')
  l=s.chomp('>').split(/[\r\n]/)
  desc=l.shift
  s=l.join("")
  r=dishuffle(s)
  puts ">#{desc} (di-shuffled)",r
end

