#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"
#include <chrono>
#include <gmock/gmock.h>

template <typename Func> auto measureTime(Func &&func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
      .count();
}

TEST(RabinKarp, SameText) {
  std::string fileName = DATASETS_PATH "1.1.txt";

  auto duration = measureTime([&]() {
    std::vector<int> result = dnaAnalyzer::dnaSearch(fileName, fileName);
    int resultSize = result.size();

    EXPECT_EQ(resultSize, 1);
    EXPECT_EQ(result[0], 0);
  });
  EXPECT_LE(duration, 1000000);
}

TEST(RabinKarp, SameTextLong) {
  std::string fileName = DATASETS_PATH "HIV-1_AF033819.3.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName, fileName);
  int resultSize = result.size();

  EXPECT_EQ(resultSize, 1);
  EXPECT_EQ(result[0], 0);
}

TEST(RabinKarp, MonoText) {
  std::string fileName1 = DATASETS_PATH "2.1.txt";
  std::string fileName2 = DATASETS_PATH "2.2.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);
  int resultSize = result.size();

  EXPECT_EQ(resultSize, 4);
  for (int i = 0; i < resultSize; i++) {
    EXPECT_EQ(result[i], i);
  }
}

TEST(RabinKarp, TextFromSubject) {
  std::string fileName1 = DATASETS_PATH "HIV-1_AF033819.3.txt";
  std::string fileName2 = DATASETS_PATH "HIV-1_substr.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);
  int resultSize = result.size();

  EXPECT_EQ(resultSize, 2);
  EXPECT_EQ(result[0], 65);
  EXPECT_EQ(result[1], 9150);
}

TEST(RabinKarp, NoMatches) {
  std::string fileName1 = DATASETS_PATH "1.2.txt";
  std::string fileName2 = DATASETS_PATH "2.2.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);

  EXPECT_TRUE(result.empty());
}

TEST(RabinKarp, NoMatchesTextLong) {
  std::string fileName1 = DATASETS_PATH "HIV-1_AF033819.3.txt";
  std::string fileName2 = DATASETS_PATH "3.2.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);

  EXPECT_TRUE(result.empty());
}

TEST(RabinKarp, EmptyPattern) {
  std::string fileName1 = DATASETS_PATH "1.1.txt";
  std::string fileName2 = DATASETS_PATH "empty.txt";

  EXPECT_THROW(dnaAnalyzer::dnaSearch(fileName1, fileName2),
               std::invalid_argument);
}

TEST(RabinKarp, EmptyText) {
  std::string fileName1 = DATASETS_PATH "empty.txt";
  std::string fileName2 = DATASETS_PATH "1.1.txt";

  EXPECT_THROW(dnaAnalyzer::dnaSearch(fileName1, fileName2),
               std::invalid_argument);
}

TEST(RabinKarp, PatternLongerThanText) {
  std::string fileName1 = DATASETS_PATH "1.2.txt";
  std::string fileName2 = DATASETS_PATH "1.1.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);

  EXPECT_TRUE(result.empty());
}

TEST(RabinKarp, WrongFile) {
  std::string fileName1 = DATASETS_PATH "nw_22.txt";
  std::string fileName2 = DATASETS_PATH "nw_2.txt";
  EXPECT_THROW(dnaAnalyzer::dnaSearch(fileName1, fileName2),
               std::invalid_argument);
}

TEST(RabinKarp, MultipleSmallPatternMatches) {
  std::string fileName1 = DATASETS_PATH "3.1.txt";
  std::string fileName2 = DATASETS_PATH "3.2.txt";
  std::vector<int> result = dnaAnalyzer::dnaSearch(fileName1, fileName2);
  int resultSize = result.size();

  EXPECT_EQ(resultSize, 6);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 3);
  EXPECT_EQ(result[2], 6);
  EXPECT_EQ(result[3], 9);
  EXPECT_EQ(result[4], 12);
  EXPECT_EQ(result[5], 15);
}
