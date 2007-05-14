#include "fasta.h"
#include <string>
#include <sstream>
#include <algorithm>    // for transform

using namespace std;

const char  HEADER = '>';
const char  TERM   = '*';
const char* BLANK  = " \t\n\r";

const int BUFSIZE = 4096;


////////////////////////////////////////////////////////////////////////////////
// Static Methods

#if 0
static string& eat(string& str, const char* blank)
{
    assert(blank);

    for (;;) {
        int pos = str.find_first_of(blank);
        if (pos == -1)
            break;
        str.erase(pos, 1);
    }
    return str;
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Constructor & Destructor

Fasta::Fasta(istream& istrm)
    : m_istrm(istrm)
{
    std::string buf;
    std::getline(m_istrm, buf, HEADER);
}

Fasta::~Fasta()
{
}


////////////////////////////////////////////////////////////////////////////////
// Interfaces

const char* Fasta::GetNextSeq()
{
    std::string buf;
    stringstream ss;

    std::getline(m_istrm, buf, HEADER);
    if (!m_istrm) {
        return "";
    }
    ss << buf;
    while (ss.get()==' ');
    ss.unget();
    std::getline(ss, m_name);
    m_seq.erase();
    while (!ss.eof()) {
        string line;
        ss >> line;
        m_seq += line;
    }
    if (*m_seq.rbegin() == TERM) {
        m_seq.erase(m_seq.length() - 1, 1);
    }

    return m_seq.c_str();
}

void Fasta::ToUpper()
{
    transform(m_seq.begin(), m_seq.end(), m_seq.begin(),
              (int (*)(int))toupper);
}

void Fasta::ToLower()
{
    transform(m_seq.begin(), m_seq.end(), m_seq.begin(),
              (int (*)(int))tolower);
}
