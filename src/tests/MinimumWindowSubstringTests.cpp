#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"
#include <gmock/gmock.h>

TEST(MinimumWindowSubstring, StringsFromSubject) {
  std::string fileName = DATASETS_PATH "minWin_1.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "GACACCCACCATACAT");
}

TEST(MinimumWindowSubstring, PatternOneSymbol) {
  std::string fileName = DATASETS_PATH "minWin_2.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "T");
}

TEST(MinimumWindowSubstring, PatternTwoMinSimbols) {
  std::string fileName = DATASETS_PATH "minWin_3.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "TCT");
}

TEST(MinimumWindowSubstring, NoResult1) {
  std::string fileName = DATASETS_PATH "minWin_4.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "");
}

TEST(MinimumWindowSubstring, Empty) {
  std::string fileName = DATASETS_PATH "empty.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "");
}

TEST(MinimumWindowSubstring, AllString) {
  std::string fileName = DATASETS_PATH "kSim_6.txt";

  EXPECT_EQ(dnaAnalyzer::minimumWindowSubstring(fileName), "ACGCGTCA");
}