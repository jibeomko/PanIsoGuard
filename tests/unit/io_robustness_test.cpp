#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>
#include <functional>
#include <string>

#include "panisoguard/bed12.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/pangenome.hpp"
#include "panisoguard/sj_tab.hpp"

using namespace panisoguard;

namespace {

// Write `content` to a unique temp file in the test working dir; removed on scope exit.
struct TempFile {
  std::string path;
  explicit TempFile(const std::string& content) {
    static int counter = 0;
    path = "./.__pig_malformed_" + std::to_string(counter++) + ".tmp";
    std::ofstream(path) << content;
  }
  ~TempFile() { std::remove(path.c_str()); }
};

// True iff `fn` throws a std::exception whose message contains every `needle`.
bool throws_containing(const std::function<void()>& fn, std::initializer_list<const char*> needles) {
  try {
    fn();
  } catch (const std::exception& e) {
    const std::string msg = e.what();
    for (const char* n : needles) {
      if (msg.find(n) == std::string::npos) return false;
    }
    return true;
  }
  return false;
}

}  // namespace

TEST_CASE("BED12: malformed numeric fields throw with file/line/field context", "[io][robustness]") {
  // chromStart not an integer
  TempFile bad_start("chr1\tXX\t500\tn\t0\t+\t100\t500\t0\t2\t100,100,\t0,400,\n");
  CHECK(throws_containing([&] { read_bed12(bad_start.path); }, {"BED12", "line 1", "chromStart"}));

  // a blockSizes token is garbage (previously a bare std::stoll terminate)
  TempFile bad_blocks("chr1\t100\t500\tn\t0\t+\t100\t500\t0\t2\t100,Q9,\t0,400,\n");
  CHECK(throws_containing([&] { read_bed12(bad_blocks.path); }, {"BED12", "blockSizes"}));
}

TEST_CASE("GTF: malformed exon coordinate throws with line context", "[io][robustness]") {
  TempFile bad("chr1\tsrc\texon\t100\tABC\t.\t+\t.\ttranscript_id \"t1\";\n");
  CHECK(throws_containing([&] { read_gtf_transcripts(bad.path); }, {"GTF", "line 1", "exon"}));
}

TEST_CASE("SJ.tab / pangenome: inverted intron (start > end) is rejected", "[io][robustness]") {
  TempFile sj("chr1\t500\t200\t1\t1\t1\t10\t0\t30\n");  // start 500 > end 200
  CHECK(throws_containing([&] { read_sj_tab(sj.path); }, {"SJ.tab", "intron_start", "intron_end"}));

  TempFile pg("chr1\t500\t200\t+\t2\n");  // start 500 > end 200
  CHECK(throws_containing([&] { read_pangenome_junctions(pg.path); },
                          {"pangenome", "intron_start", "intron_end"}));
}
