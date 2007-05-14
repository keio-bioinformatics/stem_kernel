/*!
 *  @file       fasta.h
 *  @version    2.0.0
 *  @since      2003/09/05
 *  @date       2003/12/07
 *  @author     Hiroshi Matsui

 *  @brief      FASTA parser
 */

#ifndef ___FASTA_H___
#define ___FASTA_H___

#include <istream>

class Fasta {
public:
    Fasta(std::istream& istrm);
    ~Fasta();

    // interfaces
    const char* GetNextSeq();

    const std::string& Seq() const { return m_seq; }
    const std::string& Name() const { return m_name; }

    void ToUpper();
    void ToLower();

    bool IsEOF() const { return m_istrm.eof(); }
    bool IsFail() const { return m_istrm.fail(); }
    bool operator!() const { return !m_istrm; }

protected:
    std::istream& m_istrm;

    std::string m_seq;
    std::string m_name;
};

#endif//___FASTA_H___
