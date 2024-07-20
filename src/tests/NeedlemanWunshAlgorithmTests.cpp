#include <gmock/gmock.h>
// #include "dnaAnalyzer.h"
#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"

TEST(NeedlemanWunsh, StringsFromSubjectOptimalScore) {
  std::string fileName = DATASETS_PATH "nw_1.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);

  EXPECT_EQ(optimalScore, 10);
}

TEST(NeedlemanWunsh, SameStringsOptimalScore) {
  std::string fileName = DATASETS_PATH "nw_2.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);

  EXPECT_EQ(optimalScore, 4);
}

TEST(NeedlemanWunsh, WrongFile) {
  std::string fileName = DATASETS_PATH "nw_22.txt";
  EXPECT_THROW(dnaAnalyzer::sequenceAlignmentOptimalScore(fileName),
               std::invalid_argument);
}

TEST(NeedlemanWunsh, EmptyStringsOptimalScore) {
  std::string fileName = DATASETS_PATH "nw_3_empty.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);

  EXPECT_EQ(optimalScore, 0);
}

TEST(NeedlemanWunsh, OnlyFirstStringOptimalScore) {
  std::string fileName = DATASETS_PATH "nw_4.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);

  EXPECT_EQ(optimalScore, -8);
}

TEST(NeedlemanWunsh, OnlySecondStringOptimalScore) {
  std::string fileName = DATASETS_PATH "nw_4.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);

  EXPECT_EQ(optimalScore, -8);
}

TEST(NeedlemanWunsh, StringsFromSubjectRecoveringOptimalSequenceAlignment) {
  std::string fileName = DATASETS_PATH "nw_1.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);
  std::string OptimalSequenceAlignment = "GGGCGACACTCCACCATAGA-\n"
                                         " |||||||| |||||||| | \n"
                                         "-GGCGACAC-CCACCATACAT";

  EXPECT_EQ(optimalScore, 10);
  EXPECT_EQ(OptimalSequenceAlignment,
            dnaAnalyzer::recoveringOptimalSequenceAlignment(fileName));
}

TEST(NeedlemanWunsh, SameStringsRecoveringOptimalSequenceAlignment) {
  std::string fileName = DATASETS_PATH "nw_2.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);
  std::string OptimalSequenceAlignment = "ACGT\n"
                                         "||||\n"
                                         "ACGT";

  EXPECT_EQ(optimalScore, 4);
  EXPECT_EQ(OptimalSequenceAlignment,
            dnaAnalyzer::recoveringOptimalSequenceAlignment(fileName));
}

TEST(NeedlemanWunsh, EmptyStringsRecoveringOptimalSequenceAlignment) {
  std::string fileName = DATASETS_PATH "nw_3_empty.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);
  std::string OptimalSequenceAlignment = "\n\n";

  EXPECT_EQ(optimalScore, 0);
  EXPECT_EQ(OptimalSequenceAlignment,
            dnaAnalyzer::recoveringOptimalSequenceAlignment(fileName));
}

TEST(NeedlemanWunsh, OnlyOneStringRecoveringOptimalSequenceAlignment) {
  std::string fileName = DATASETS_PATH "nw_4.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);
  std::string OptimalSequenceAlignment = "ACGT\n"
                                         "    \n"
                                         "----";

  EXPECT_EQ(optimalScore, -8);
  EXPECT_EQ(OptimalSequenceAlignment,
            dnaAnalyzer::recoveringOptimalSequenceAlignment(fileName));
}

TEST(NeedlemanWunsh, ShortStringRecoveringOptimalSequenceAlignment) {
  std::string fileName = DATASETS_PATH "nw_5.txt";
  long long optimalScore = dnaAnalyzer::sequenceAlignmentOptimalScore(fileName);
  std::string OptimalSequenceAlignment = "ACGT-\n"
                                         " ||| \n"
                                         "-CGTA";

  EXPECT_EQ(optimalScore, -1);
  EXPECT_EQ(OptimalSequenceAlignment,
            dnaAnalyzer::recoveringOptimalSequenceAlignment(fileName));
}
