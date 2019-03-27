#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/memory/memory.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "fmt/format.h"
#include "leveldb/write_batch.h"
#include "roaring64map.hh"

#include "ggff.h"

namespace ggff {
void populate_subgraph(const std::string &ggff_col1, pb::SubGraph *sg) {
  // {NODE_ID}[{NODE_START}:{NODE_END}][+-],...
  const std::vector<std::string> tokens =
      absl::StrSplit(ggff_col1, absl::ByChar(','), absl::SkipEmpty());

  const auto toI64 = [](unsigned long long int val) -> std::uint64_t {
    return static_cast<std::uint64_t>(val);
  };

  for (const auto &token : tokens) {
    const std::vector<std::string> nodeTokens =
        absl::StrSplit(token, absl::ByAnyChar("[:]"), absl::SkipEmpty());

    if (nodeTokens.size() != 4) {
      std::cerr << ggff_col1 << std::endl;
    }

    sg->add_ids(toI64(std::strtoull(nodeTokens[0].c_str(), nullptr, 10)));
    auto slice = sg->add_offsets();
    slice->set_start(toI64(std::strtoull(nodeTokens[1].c_str(), nullptr, 10)));
    slice->set_end(toI64(std::strtoull(nodeTokens[2].c_str(), nullptr, 10)));
    sg->add_strands(nodeTokens[3] == "?"
                        ? pb::UNSPECIFIED_STRAND
                        : nodeTokens[3] == "-" ? pb::REVERSE_STRAND
                                               : pb::FORWARD_STRAND);
  }
}

pb::GGFF parse_ggff(const std::string &ggff_path) {
  std::ifstream fileHandle(ggff_path, std::ios_base::in);
  std::string ggffLine;

  pb::GGFF result;
  while (std::getline(fileHandle, ggffLine)) {
    std::vector<std::string> tokens =
        absl::StrSplit(ggffLine, absl::ByChar('\t'), absl::SkipEmpty());

    auto rec = result.add_records();
    populate_subgraph(tokens[0], rec->mutable_sub_graph());
    rec->set_ggff_source(tokens[1]);
    rec->set_ggff_type(tokens[2]);

    tokens.erase(tokens.begin());
    rec->set_ggff_line(absl::StrJoin(tokens.cbegin(), tokens.cend(), "\t"));
  }

  return result;
}

void write_ggff(const pb::GGFF &g, std::ostream &oss) {
  std::string result;
  for (const auto &rec : g.records()) {
    const auto &sg = rec.sub_graph();
    const auto numNodes = sg.ids_size();

    for (auto idx = 0; idx < numNodes; idx++) {
      result +=
          fmt::format("{}[{}:{}]{}", sg.ids(idx), sg.offsets(idx).start(),
                      sg.offsets(idx).end(),
                      sg.strands(idx) == pb::UNSPECIFIED_STRAND
                          ? '?'
                          : sg.strands(idx) == pb::REVERSE_STRAND ? '-' : '+');
    }

    if (!result.empty()) {
      result.pop_back();
    }

    result += fmt::format("\t{}\n", rec.ggff_line());
  }

  oss << result;
}

void write_subgraphs(const std::vector<pb::SubGraph> &graphs,
                     std::ostream &oss) {
  std::string result;
  for (const auto &sg : graphs) {
    const auto numNodes = sg.ids_size();

    for (auto idx = 0; idx < numNodes; idx++) {
      result +=
          fmt::format("{}[{}:{}]{}", sg.ids(idx), sg.offsets(idx).start(),
                      sg.offsets(idx).end(),
                      sg.strands(idx) == pb::UNSPECIFIED_STRAND
                          ? '?'
                          : sg.strands(idx) == pb::REVERSE_STRAND ? '-' : '+');
    }

    if (!result.empty()) {
      result.pop_back();
      result += '\n';
    }
  }

  oss << result;
}

std::string record_key(const pb::GGFFRecord &rec) {
  return fmt::format("{}:{}", rec.ggff_type(), rec.ggff_source());
}

inline absl::flat_hash_map<std::uint64_t, std::size_t>
build_node2offset_idx(const pb::SubGraph &sg) {
  absl::flat_hash_map<std::uint64_t, std::size_t> node2OffsetIdx;
  std::size_t idx = 0;
  for (const auto &nodeID : sg.ids()) {
    node2OffsetIdx[nodeID] = idx;
    idx++;
  }
  return node2OffsetIdx;
}

pb::Strand operator&(const pb::Strand &s1, const pb::Strand &s2) {
  return s1 == s2 ? s1 : pb::UNSPECIFIED_STRAND;
}

pb::Strand operator-(const pb::Strand &s1, const pb::Strand &s2) {
  return s1 != s2 ? s1 : pb::UNSPECIFIED_STRAND;
}

pb::Strand operator|(const pb::Strand &s1, const pb::Strand &s2) {
  return s1 == s2 ? s1 : pb::UNSPECIFIED_STRAND;
}

pb::Strand operator^(const pb::Strand &s1, const pb::Strand &s2) {
  return s1 != s2 ? s1 : pb::UNSPECIFIED_STRAND;
}

inline Roaring64Map build_slice_bitmap(const pb::NodeSlice &sl) {
  Roaring64Map result;
  for (auto idx = sl.start(); idx < sl.end(); idx++) {
    result.add(idx);
  }
  return result;
}

pb::NodeSlice operator&(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2) {
  const auto resultMap = build_slice_bitmap(sl1) & build_slice_bitmap(sl2);
  pb::NodeSlice result;
  result.set_start(resultMap.minimum());
  result.set_end(resultMap.maximum());
  return result;
}

pb::NodeSlice operator-(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2) {
  const auto resultMap = build_slice_bitmap(sl1) - build_slice_bitmap(sl2);
  pb::NodeSlice result;
  result.set_start(resultMap.minimum());
  result.set_end(resultMap.maximum());
  return result;
}

pb::NodeSlice operator|(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2) {
  const auto resultMap = build_slice_bitmap(sl1) | build_slice_bitmap(sl2);
  pb::NodeSlice result;
  result.set_start(resultMap.minimum());
  result.set_end(resultMap.maximum());
  return result;
}

pb::NodeSlice operator^(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2) {
  const auto resultMap = build_slice_bitmap(sl1) ^ build_slice_bitmap(sl2);
  pb::NodeSlice result;
  result.set_start(resultMap.minimum());
  result.set_end(resultMap.maximum());
  return result;
}

pb::SubGraph operator&(const pb::SubGraph &sg1, const pb::SubGraph &sg2) {
  // build nodeID -> offset index map for each sub graph
  const auto sg1OffsetIdx = build_node2offset_idx(sg1);
  const auto sg2OffsetIdx = build_node2offset_idx(sg2);

  const Roaring64Map sg1NodeMap(static_cast<std::size_t>(sg1.ids_size()),
                                sg1.ids().data());
  const Roaring64Map sg2NodeMap(static_cast<std::size_t>(sg2.ids_size()),
                                sg2.ids().data());

  const auto resultNodeMap = sg1NodeMap & sg2NodeMap;
  const auto resultSize = resultNodeMap.cardinality();
  auto nodeIDs = new std::uint64_t[resultSize];

  pb::SubGraph out;
  for (auto idx = 0; idx < resultSize; idx++) {
    out.add_ids(nodeIDs[idx]);

    const auto sg1Idx = static_cast<int>(sg1OffsetIdx.at(nodeIDs[idx]));
    const auto sg2Idx = static_cast<int>(sg2OffsetIdx.at(nodeIDs[idx]));

    auto slice = out.add_offsets();
    const auto sliceResult = sg1.offsets(sg1Idx) & sg2.offsets(sg2Idx);
    slice->set_start(sliceResult.start());
    slice->set_end(sliceResult.end());
    out.add_strands(sg1.strands(sg1Idx) & sg2.strands(sg2Idx));
  }

  delete[] nodeIDs;
  return out;
}

pb::SubGraph operator-(const pb::SubGraph &sg1, const pb::SubGraph &sg2) {
  // build nodeID -> offset index map for each sub graph
  const auto sg1OffsetIdx = build_node2offset_idx(sg1);
  const auto sg2OffsetIdx = build_node2offset_idx(sg2);

  const Roaring64Map sg1NodeMap(static_cast<std::size_t>(sg1.ids_size()),
                                sg1.ids().data());
  const Roaring64Map sg2NodeMap(static_cast<std::size_t>(sg2.ids_size()),
                                sg2.ids().data());

  const auto resultNodeMap = sg1NodeMap - sg2NodeMap;
  const auto resultSize = resultNodeMap.cardinality();
  auto nodeIDs = new std::uint64_t[resultSize];

  pb::SubGraph out;
  for (auto idx = 0; idx < resultSize; idx++) {
    out.add_ids(nodeIDs[idx]);

    const auto sg1Idx = static_cast<int>(sg1OffsetIdx.at(nodeIDs[idx]));
    const auto sg2Idx = static_cast<int>(sg2OffsetIdx.at(nodeIDs[idx]));

    auto slice = out.add_offsets();
    const auto sliceResult = sg1.offsets(sg1Idx) - sg2.offsets(sg2Idx);
    slice->set_start(sliceResult.start());
    slice->set_end(sliceResult.end());
    out.add_strands(sg1.strands(sg1Idx) - sg2.strands(sg2Idx));
  }

  delete[] nodeIDs;
  return out;
}

pb::SubGraph operator|(const pb::SubGraph &sg1, const pb::SubGraph &sg2) {
  // build nodeID -> offset index map for each sub graph
  const auto sg1OffsetIdx = build_node2offset_idx(sg1);
  const auto sg2OffsetIdx = build_node2offset_idx(sg2);

  const Roaring64Map sg1NodeMap(static_cast<std::size_t>(sg1.ids_size()),
                                sg1.ids().data());
  const Roaring64Map sg2NodeMap(static_cast<std::size_t>(sg2.ids_size()),
                                sg2.ids().data());

  const auto resultNodeMap = sg1NodeMap | sg2NodeMap;
  const auto resultSize = resultNodeMap.cardinality();
  auto nodeIDs = new std::uint64_t[resultSize];

  pb::SubGraph out;
  for (auto idx = 0; idx < resultSize; idx++) {
    out.add_ids(nodeIDs[idx]);

    const auto sg1Idx = static_cast<int>(sg1OffsetIdx.at(nodeIDs[idx]));
    const auto sg2Idx = static_cast<int>(sg2OffsetIdx.at(nodeIDs[idx]));

    auto slice = out.add_offsets();
    const auto sliceResult = sg1.offsets(sg1Idx) | sg2.offsets(sg2Idx);
    slice->set_start(sliceResult.start());
    slice->set_end(sliceResult.end());
    out.add_strands(sg1.strands(sg1Idx) | sg2.strands(sg2Idx));
  }

  delete[] nodeIDs;
  return out;
}

pb::SubGraph operator^(const pb::SubGraph &sg1, const pb::SubGraph &sg2) {
  // build nodeID -> offset index map for each sub graph
  const auto sg1OffsetIdx = build_node2offset_idx(sg1);
  const auto sg2OffsetIdx = build_node2offset_idx(sg2);

  const Roaring64Map sg1NodeMap(static_cast<std::size_t>(sg1.ids_size()),
                                sg1.ids().data());
  const Roaring64Map sg2NodeMap(static_cast<std::size_t>(sg2.ids_size()),
                                sg2.ids().data());

  const auto resultNodeMap = sg1NodeMap ^ sg2NodeMap;
  const auto resultSize = resultNodeMap.cardinality();
  auto nodeIDs = new std::uint64_t[resultSize];

  pb::SubGraph out;
  for (auto idx = 0; idx < resultSize; idx++) {
    out.add_ids(nodeIDs[idx]);

    const auto sg1Idx = static_cast<int>(sg1OffsetIdx.at(nodeIDs[idx]));
    const auto sg2Idx = static_cast<int>(sg2OffsetIdx.at(nodeIDs[idx]));

    auto slice = out.add_offsets();
    const auto sliceResult = sg1.offsets(sg1Idx) ^ sg2.offsets(sg2Idx);
    slice->set_start(sliceResult.start());
    slice->set_end(sliceResult.end());
    out.add_strands(sg1.strands(sg1Idx) ^ sg2.strands(sg2Idx));
  }

  delete[] nodeIDs;
  return out;
}

GGFFIndex::GGFFIndex(const std::string &ggff_index_path) {
  leveldb::Options opts;
  opts.block_size = 1 << 20; // 1 MB blocks
  opts.create_if_missing = true;

  leveldb::DB *db = nullptr;
  const auto status = leveldb::DB::Open(opts, ggff_index_path, &db);
  assert(status.ok());
  index = absl::WrapUnique(db);
}

void GGFFIndex::add_ggff(const std::string &ggff_path) {
  const auto GGFF = parse_ggff(ggff_path);

  leveldb::WriteBatch batch;
  std::string valueBuffer;
  for (const auto &rec : GGFF.records()) {
    rec.SerializeToString(&valueBuffer);
    batch.Put(record_key(rec), valueBuffer);
  }

  const auto status = index->Write(leveldb::WriteOptions(), &batch);
  assert(status.ok());
}

void GGFFIndex::iterate(GGFFIndex::HandleRecord &handle_func,
                        const std::string &record_key_prefix) {
  auto itr = absl::WrapUnique(index->NewIterator(leveldb::ReadOptions()));
  record_key_prefix.empty() ? itr->SeekToFirst() : itr->Seek(record_key_prefix);
  pb::GGFFRecord rec;
  for (; itr->Valid(); itr->Next()) {
    assert(rec.ParseFromString(itr->value().ToString()));
    handle_func(rec);
  }
}

std::vector<pb::SubGraph>
GGFFIndex::intersect_of(const pb::SubGraph &g,
                        const std::string &record_key_prefix) {
  std::vector<pb::SubGraph> results;
  HandleRecord handler = [&g, &results](const pb::GGFFRecord &current) {
    const auto result = g & current.sub_graph();
    if (!result.ids().empty()) {
      results.push_back(result);
    }
  };

  iterate(handler, record_key_prefix);
  return results;
}

std::vector<pb::SubGraph>
GGFFIndex::difference_of(const pb::SubGraph &g,
                         const std::string &record_key_prefix) {
  std::vector<pb::SubGraph> results;
  HandleRecord handler = [&g, &results](const pb::GGFFRecord &current) {
    const auto result = g - current.sub_graph();
    if (!result.ids().empty()) {
      results.push_back(result);
    }
  };

  iterate(handler, record_key_prefix);
  return results;
}

std::vector<pb::SubGraph>
GGFFIndex::union_of(const pb::SubGraph &g,
                    const std::string &record_key_prefix) {
  std::vector<pb::SubGraph> results;
  HandleRecord handler = [&g, &results](const pb::GGFFRecord &current) {
    const auto result = g | current.sub_graph();
    if (!result.ids().empty()) {
      results.push_back(result);
    }
  };

  iterate(handler, record_key_prefix);
  return results;
}

std::vector<pb::SubGraph>
GGFFIndex::xor_of(const pb::SubGraph &g, const std::string &record_key_prefix) {
  std::vector<pb::SubGraph> results;
  HandleRecord handler = [&g, &results](const pb::GGFFRecord &current) {
    const auto result = g ^ current.sub_graph();
    if (!result.ids().empty()) {
      results.push_back(result);
    }
  };

  iterate(handler, record_key_prefix);
  return results;
}
} // namespace ggff