#!/usr/bin/env ruby

def roc(x)
  pos=[]
  neg=[]
  x.each do |c,v|
    if c>=0
      pos.push(v)
    else
      neg.push(v)
    end
  end
  pos.sort!{|a,b| b<=>a}
  neg.sort!{|a,b| b<=>a}

  cur=[]
  tp=fp=0
  i=j=0
  while i<pos.size && j<neg.size
    cur.push([fp/neg.size.to_f, tp/pos.size.to_f])
    if pos[i]>neg[j]
      tp+=1; i+=1;
    elsif pos[i]<neg[j]
      fp+=1; j+=1;
    else
      fp+=1; j+=1;
      tp+=1; i+=1;
    end
  end
  cur.push([fp/neg.size.to_f, tp/pos.size.to_f])
  cur.push([1.0, 1.0])
  
  s=0.0
  (1..cur.size-1).each do|i|
    s += (cur[i-1][1]+cur[i][1])*(cur[i][0]-cur[i-1][0])/2
  end
  [s,cur]
end

def accspsn(x, th=0.0)
  tp=fp=tn=fn=0
  x.each do |c,v|
    if c>=0
      if v>=th
        tp+=1
      else
        fn+=1
      end
    else
      if v>=th
        fp+=1
      else
        tn+=1
      end
    end
  end
  [(tp+tn).to_f/(tp+tn+fp+fn),
   tn.to_f/(tn+fp),
   tp.to_f/(tp+fn)]
end

if $0=~/roc.rb/
  x=[]
  while l=gets
    c,v=l.chomp.split
    x.push([c.to_i,v.to_f])
  end
  s,r=roc(x)
  acc,sp,sn=accspsn(x)
  puts "acc=#{acc*100}, sp=#{sp*100}, sn=#{sn*100}"
  puts "ROC score = #{s}"
#   r.each{|x|
#     puts x.join(" ")
#   }
end
