%module onposix
%{
#include "Assert.hpp"
#include "AbstractDescriptorReader.hpp"
#include "AbstractThread.hpp"
#include "PosixCondition.hpp"
#include "PosixDescriptor.hpp"
#include "PosixMutex.hpp"
#include "PosixPrioritySharedQueue.hpp"
#include "PosixSharedQueue.hpp"
#include "Buffer.hpp"
#include "Process.hpp"
#include "SimpleThread.hpp"
#include "StreamSocketClientDescriptor.hpp"
#include "StreamSocketServerDescriptor.hpp"
#include "StreamSocketServer.hpp"
#include "Time.hpp"
#include "DescriptorsMonitor.hpp"
#include "DgramSocketClientDescriptor.hpp"
#include "DgramSocketServerDescriptor.hpp"
#include "FifoDescriptor.hpp"
#include "FileDescriptor.hpp"
#include "Logger.hpp"
#include "Pipe.hpp"

using namespace onposix;
%}

%include "PosixCondition.hpp"
%include "PosixDescriptor.hpp"
%include "PosixMutex.hpp"
%include "PosixPrioritySharedQueue.hpp"
%include "PosixSharedQueue.hpp"
%include "Buffer.hpp"
%include "Process.hpp"
%include "SimpleThread.hpp"
%include "StreamSocketClientDescriptor.hpp"
%include "StreamSocketServerDescriptor.hpp"
%include "StreamSocketServer.hpp"
%include "Time.hpp"
%include "DescriptorsMonitor.hpp"
%include "DgramSocketClientDescriptor.hpp"
%include "DgramSocketServerDescriptor.hpp"
%include "FifoDescriptor.hpp"
%include "FileDescriptor.hpp"
%include "Logger.hpp"
%include "Pipe.hpp"

%extend onposix::Buffer {

    char __getitem(size_t i) { return (*$self)[i]; }
    void __setitem(size_t i, char c) { (*$self)[i] = c; }
}