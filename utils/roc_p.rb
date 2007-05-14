#!/usr/bin/env ruby

require 'roc'

pos_label=1

ans=[]
IO.foreach(ARGV.shift) do |l|
  a,=l.chomp.split
  ans.push a.to_i
end

l=gets
l=l.chomp.split
l.map!{|v| v.to_i}
pos=l.index(pos_label)
raise if pos<0

tp=fp=tn=fn=0
v=[]
n=0
while l=gets
  l=l.chomp.split
  v.push l[pos].to_f

  if l[0].to_i==ans[n]
    if ans[n]>=0; tp+=1; else; tn+=1; end
  else
    if ans[n]>=0; fn+=1; else; fp+=1; end
  end
  n+=1
end

raise unless ans.size==v.size
x=[]
ans.each_index{|i| x.push [ans[i],v[i]]}
s,r=roc(x)
acc=(tp+tn)/n.to_f
sp=tn/(tn+fp).to_f
sn=tp/(tp+fn).to_f
puts "acc=#{acc*100}, sp=#{sp*100}, sn=#{sn*100}"
puts "ROC score = #{s}"
