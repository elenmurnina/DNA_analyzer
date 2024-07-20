#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"
#include <gmock/gmock.h>

TEST(MatchingRegularExpression, StringsFromSubject) {
  std::string fileName = DATASETS_PATH "regExp_1.txt";

  EXPECT_EQ(dnaAnalyzer::matchingRegularExpressions(fileName), true);
}

TEST(MatchingRegularExpression, NonMatchingPatternFalse) {
  std::string fileName = DATASETS_PATH "regExp_2.txt";

  EXPECT_EQ(dnaAnalyzer::matchingRegularExpressions(fileName), false);
}

TEST(MatchingRegularExpression, Empty) {
  std::string fileName = DATASETS_PATH "empty.txt";

  EXPECT_EQ(dnaAnalyzer::matchingRegularExpressions(fileName), true);
}

TEST(MatchingRegularExpression, MatchingPatternTrue1) {
  std::string fileName = DATASETS_PATH "regExp_3.txt";

  EXPECT_EQ(dnaAnalyzer::matchingRegularExpressions(fileName), true);
}

TEST(MatchingRegularExpression, InvalidCharacterInSequence) {
  std::string fileName = DATASETS_PATH "regExp_4.txt";

  EXPECT_THROW(dnaAnalyzer::matchingRegularExpressions(fileName),
               std::invalid_argument);
}

TEST(MatchingRegularExpression, InvalidCharacterInPattern) {
  std::string fileName = DATASETS_PATH "regExp_5.txt";

  EXPECT_THROW(dnaAnalyzer::matchingRegularExpressions(fileName),
               std::invalid_argument);
}

TEST(MatchingRegularExpression, InvalidCharacterInPattern1) {
  std::string fileName = DATASETS_PATH "regExp_6.txt";

  EXPECT_THROW(dnaAnalyzer::matchingRegularExpressions(fileName),
               std::invalid_argument);
}

TEST(MatchingRegularExpression, MatchingPatternTrue2) {
  std::string fileName = DATASETS_PATH "regExp_7.txt";

  EXPECT_EQ(dnaAnalyzer::matchingRegularExpressions(fileName), true);
}