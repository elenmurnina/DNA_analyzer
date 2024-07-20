#include "../dnaAnalyzerImpl/includes/dnaAnalyzer.h"
#include <gmock/gmock.h>

TEST(dnaAnalyzerTest, ConstructorTest) { dnaAnalyzer analyzer; }

TEST(dnaAnalyzerTest, DestructorTest) {
  dnaAnalyzer *analyzer = new dnaAnalyzer();
  delete analyzer;
}
