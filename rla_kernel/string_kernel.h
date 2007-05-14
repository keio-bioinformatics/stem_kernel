
// $Id: string_kernel.h 64 2006-08-02 02:17:48Z satoken $

#ifndef __INC_STRING_KERNEL_H__
#define __INC_STRING_KERNEL_H__

#include <string>
#include <vector>

template < class ValueType >
class StringKernel
{
 public:
  typedef ValueType value_type;
  
  class Seq{
  public:
    std::string seq;
    std::vector<value_type> p_r;
    std::vector<value_type> p_l;
    //std::vector<value_type> p_un;
    //std::vector<uint> seq;
    Seq() : seq(), p_r(), p_l(){}
    Seq(const std::string& s, 
	const std::vector<value_type>& r,
	const std::vector<value_type>& l)
      : seq(s), p_r(r), p_l(l){}
    //Seq() : seq(), p_r(), p_l(), p_un(){}
    //Seq(const std::vector<uint>& s, std::vector<value_type> r, std::vector<value_type> l,std::vector<value_type> un): seq(s), p_r(r), p_l(l), p_un(un){}
    //Seq(const Seq& x) : seq(x.seq), p_r(x.p_r), p_l(x.p_l), p_un(x.p_un) {}
    Seq(const Seq& x) : seq(x.seq), p_r(x.p_r), p_l(x.p_l){}
    uint size()const{ return seq.size();}
  };

 private:
  value_type gap_;
  value_type ext_;
  value_type alpha_;
  value_type beta_;
  
public:
  StringKernel(value_type gap=1,value_type ext=1,
	       value_type alpha=1,value_type beta=1) 
    : gap_(gap), ext_(ext), alpha_(alpha), beta_(beta) {}
  value_type operator()(const Seq& xx, const Seq& yy) const;
  value_type operator()(const std::string& x, const std::string& y) const;
  value_type score(const Seq& xx, const Seq& yy, uint i, uint j) const;
  int char2rna(char c) const;
};

#endif
