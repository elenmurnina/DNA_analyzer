#ifndef DNAANALYZER_H
#define DNAANALYZER_H

#include <string>
#include <vector>

class dnaAnalyzer {

public:
  dnaAnalyzer();

  static std::vector<int> dnaSearch(std::string &firstFilename,
                                    std::string &secondFilename);
  static long long sequenceAlignmentOptimalScore(std::string &filename);
  static std::string recoveringOptimalSequenceAlignment(std::string &filename);
  static bool matchingRegularExpressions(std::string &filename);
  static int kSimilarityOfTwoSequences(std::string &filename);
  static std::string minimumWindowSubstring(std::string &filename);

  ~dnaAnalyzer();
};

#endif