CXX = g++
CXXFLAGS = -O3 -std=c++11
HEADERS = alignGraphs.hpp alignmentExtender.hpp blastDistance.hpp computeSpectralSignatures.hpp config.hpp dalpha.hpp graph.hpp gzReader.hpp gzWriter.hpp hungarian.hpp laplacian.hpp progressBar.hpp readFromGexf.hpp spectralToDistance.hpp swapEndian.hpp
BOOSTDIR ?= /usr/local/include
BOOSTLIB ?= /usr/local/lib
LFLAGS = -lboost_iostreams -lboost_system -lboost_thread -lpthread

GHOST: GHOST.cpp $(HEADERS)
	$(CXX) GHOST.cpp $(CXXFLAGS) -I $(BOOSTDIR) -L $(BOOSTLIB) $(LFLAGS) -o GHOST

likenew:
	\rm -f *.sig.gz *.sdf *.af GHOST

clean:
	\rm -f GHOST

tar:
	tar cfv *.sig.gz *.sdf *.af
