#!/usr/bin/env ruby

def arrival?(x,f,edge)
  return true if x==f
  arrived=Hash.new
  until arrived[x]
    return true if edge[x][0]==f
    arrived[x]=true
    x=edge[x][0]
  end
  return false
end

def check(f,ledge)
  ledge.each{|x,y|
    return false unless arrival?(x,f,ledge)
  }
  return true
end

def make_last_graph(f,edge)
  ledge={}
  edge.each{|x,l|
    next if x==f
    i=rand(l.size)
    ledge[x]=l[i]
    l.delete_at(i)
  }
  ledge
end

def merge_last_graph(edge,ledge)
  ledge.each{|x,y| edge[x].push y }
end

def dishuffle(s, n=100)
  if s.kind_of?(String)
    dishuffle_array(s.split(''))[0].join('')
  else
    dishuffle_array(s)
  end
end

def dishuffle_array(s, n=100)
  # make Eulerian edge ordering
  edge={}
  (1..s.size-1).each{|i|
    edge[s[i-1]]=[] if edge[s[i-1]].nil?
    edge[s[i-1]].push([s[i],i])
  }
  
  # make a last-egdge graph randomly
  ledge = make_last_graph(s[-1], edge)
  until check(s[-1], ledge)
    merge_last_graph(edge, ledge)
    ledge = make_last_graph(s[-1], edge)
  end

  # randomize the Eulerian edge ordering
  edge.each{|x,l|
    (l.size*n).times do
      a=rand(l.size)
      b=rand(l.size)
      l[a],l[b]=l[b],l[a]
    end
  }

  # make string from the Eulerian edge ordering
  merge_last_graph(edge, ledge)
  ret=[[s[0],0]]
  begin
    #puts ret[-1]
    while !edge[ret[-1][0]].nil? && edge[ret[-1][0]].size>0
      #p edge
      #puts "#{ret[-1]} -> #{edge[ret[-1]]}"
      ret.push(edge[ret[-1][0]].shift)
    end
  rescue 
    raise
  end
  [ret.map{|v| v[0]}, ret.map{|v| v[1]}]
end
