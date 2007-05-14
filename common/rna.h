// $Id$

#ifndef __INC_RNA_H__
#define __INC_RNA_H__

#include <vector>
#include <deque>
#include <string>

enum {
  RNA_A = 0,
  RNA_C = 1,
  RNA_G = 2,
  RNA_T = 3,
  RNA_U = 3,
  RNA_GAP = 4,
  N_RNA = 4
};
typedef unsigned char rna_t;

typedef std::vector<rna_t> RNASequence;

template <class TP>
struct RNASymbol {
  enum {
    A = RNA_A,
    C = RNA_C,
    G = RNA_G,
    T = RNA_T,
    U = RNA_U,
    GAP = RNA_GAP,
  };
};

template <>
struct RNASymbol<char> {
  enum {
    A = 'a',
    C = 'c',
    G = 'g',
    U = 'u',
    T = 't',
    GAP = '-',
  };
};

rna_t char2rna(char rna);
char rna2char(const rna_t& rna);
void char2rna(RNASequence& out, const std::string& in);
void char2rna(std::string& out, const std::string& in);
void rna2char(std::string& out, const RNASequence& in);
void rna2char(std::string& out, const std::string& in);

template <class T>
inline
rna_t
index(const T& x)
{
  return char2rna(x);
}

template <>
inline
rna_t
index(const rna_t& x)
{
  return x;
}

// for multiple alignments

template < class V >
class Column {
public:
  typedef V value_type;
  
  Column()
    : val_(), cnt_(N_RNA, 0)
  {
  }

  Column(const Column& col)
    : val_(col.val_), cnt_(col.cnt_)
  {
  }

  Column(const std::deque<V>& val);

  Column& operator=(const Column& col)
  {
    if (this != &col) {
      val_ = col.val_;
      cnt_ = col.cnt_;
    }
    return *this;
  }

  uint n_seqs() const { return val_.size(); }
  const V& operator[](uint i) const { return val_[i]; }
  void push_back(const V& v);
  uint cnt(uint x) const { return cnt_[x]; }

private:
  std::deque<V> val_;
  std::vector<unsigned char> cnt_;
};

template < class Seq >
class MASequence : public std::vector< Column<typename Seq::value_type>  >
{
public:
  typedef Seq SingleSeq;

  uint n_seqs() const
  {
    return this->size()==0 ? 0 : (*this)[0].n_seqs();
  }

  MASequence<Seq>& add_seq(const Seq& seq);
  Seq get_seq(uint i) const;
};


template < class Seq >
Column<typename Seq::value_type>
gap_symbol(const MASequence<Seq>& ma_seq);

template <class Seq>
Seq erase_gap(const Seq& seq);

template <class Seq>
void
erase_gap(const Seq& seq, const std::string& str,
	  Seq& r_seq, std::string& r_str);

#endif	// __INC_RNA_H__

// Local Variables:
// mode: C++
// End:
