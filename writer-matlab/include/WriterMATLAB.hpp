#ifndef WRITERMATLAB_H
#define WRITERMATLAB_H

#include <writer-matlab_EXPORTS.hpp>
#include <Writer.hpp>

namespace mbsolve {

class WRITER_MATLAB_EXPORT WriterMATLAB : public IWriter
{
public:
    std::string getExtension() const;

    void write(const std::string& file, const std::vector<Result *>& results,
	       const Device& device, const Scenario& scenario) const;

};

}

#endif
