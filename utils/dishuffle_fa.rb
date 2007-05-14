#!/usr/bin/env ruby

require 'dishuffle'

gets('>')
while s=gets('>')
  l=s.split("\n")
  desc=l.shift
  s=l.join("").delete(">")
  r=dishuffle(s)
  puts ">#{desc} (di-shuffled)",r
end

