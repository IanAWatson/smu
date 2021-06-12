#include "smu.h"

namespace smu {

SmuResults::SmuResults() {
  _molecules_examined = 0;
  _duplicates_discarded = 0;
  _invalid_valence_rejected = 0;
  _separator = ',';
}

SmuResults::~SmuResults() {
}

size_t
SmuResults::size() const {
  return _smiles_to_id.size();
}

void
SmuResults::Report(std::ostream& output) const {
  output << "Examined " << _molecules_examined << " molecules\n";
  output << _duplicates_discarded << " duplicates discarded\n";
  output << _invalid_valence_rejected << " discarded for invalid valence\n";
}

bool
SmuResults::Add(std::unique_ptr<Molecule>& mol,
                bool non_aromatic_unique_smiles) {
  _molecules_examined++;
  if (_molecules_examined % 10000 == 0) {
    std::cerr << "examined " << _molecules_examined << " molecules\n";
  }

  if (! mol->valence_ok()) {
    std::cerr << " invalid valence " << mol->smiles() << "\n";
    _invalid_valence_rejected++;
    return 0;
  }

  const IWString& usmi = non_aromatic_unique_smiles ?
        mol->non_aromatic_unique_smiles() : mol->unique_smiles();

  if (_smiles_to_id.contains(usmi)) {
    _duplicates_discarded++;
    return false;
  }
  size_t id = _smiles_to_id.size() + 1;
  _smiles_to_id.emplace(usmi, id);

  return true;
}


// If the `m` passes quality tests, write to `output`.
// Returns true if written.
bool
SmuResults::_maybe_write_molecule(Molecule& m,
                   const IWString& smiles,
                   const IWString& id,
                   std::ostream& output) const {
  if (! m.valence_ok()) {
    return false;
  }

  if (! OkCharges(m)) {
    return false;
  }

  output << smiles << _separator << id << '\n';
  return true;
}

int
SmuResults::Write(std::ostream& output) const {
  int written = 0;
  for (const auto& [smiles, id] : _smiles_to_id) {
    Molecule m;
    if (! m.build_from_smiles(smiles)) {
      cerr << "Cannot build from " << smiles << "\n";
      continue;
    }
    if (_maybe_write_molecule(m, smiles, id, output)) {
      written++;
    }
  }
  return written;
}

std::vector<IWString>
SmuResults::CurrentSmiles() const {
  std::vector<IWString> result;  // to be returned.
  result.reserve(_smiles_to_id.size());
  for (auto& [smiles, mol] : _smiles_to_id) {
    result.emplace_back(smiles);
  }

  return result;
}

}  // namespace smu
