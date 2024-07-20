#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"
#include <gmock/gmock.h>

TEST(kSimilarityOfTwoSequences, StringsFromSubject) {
  std::string fileName = DATASETS_PATH "kSim_1.txt";

  EXPECT_EQ(dnaAnalyzer::kSimilarityOfTwoSequences(fileName), 3);
}

TEST(kSimilarityOfTwoSequences, StringsNotAnagrams) {
  std::string fileName = DATASETS_PATH "kSim_2.txt";

  EXPECT_THROW(dnaAnalyzer::kSimilarityOfTwoSequences(fileName),
               std::invalid_argument);
}

TEST(kSimilarityOfTwoSequences, Empty) {
  std::string fileName = DATASETS_PATH "empty.txt";

  EXPECT_EQ(dnaAnalyzer::kSimilarityOfTwoSequences(fileName), 0);
}

TEST(kSimilarityOfTwoSequences, IdenticaStrings) {
  std::string fileName = DATASETS_PATH "kSim_3.txt";

  EXPECT_EQ(dnaAnalyzer::kSimilarityOfTwoSequences(fileName), 0);
}

TEST(kSimilarityOfTwoSequences, Strings1) {
  std::string fileName = DATASETS_PATH "kSim_4.txt";

  EXPECT_EQ(dnaAnalyzer::kSimilarityOfTwoSequences(fileName), 3);
}

TEST(kSimilarityOfTwoSequences, Strings2) {
  std::string fileName = DATASETS_PATH "kSim_6.txt";

  EXPECT_EQ(dnaAnalyzer::kSimilarityOfTwoSequences(fileName), 5);
}

TEST(kSimilarityOfTwoSequences, DiffStrings) {
  std::string fileName = DATASETS_PATH "kSim_5.txt";

  EXPECT_THROW(dnaAnalyzer::kSimilarityOfTwoSequences(fileName),
               std::invalid_argument);
}