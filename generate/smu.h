#ifndef MOLECULE_TOOLS_SMU_H_
#define MOLECULE_TOOLS_SMU_H_

#include <memory>
#include <vector>

#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"

namespace smu {

struct SMU_Params {
  resizable_array_p<Atom> to_add;
  int max_atoms;
  bool non_aromatic_unique_smiles;
  bool oneg_only_added_to_npos;
};

struct MoleculeData {
  // For each atom, the max number of connections.
  std::vector<int> max_conn;
  // For each atom, the max number of bonds.
  std::vector<int> max_bonds;
};

// Results can be stored in a large or small memory form.
class SmuResults {
  public:
    SmuResults();
    ~SmuResults();

    bool Add(std::unique_ptr<Molecule>& mol, bool non_aromatic_unique_smiles);

    void set_output_separator(char s) {
      _separator = s;
    }

    // Report statistics about molecules examined.
    void Report(std::ostream& output) const;

    // Write smiles file to `output`.
    int Write(std::ostream& output) const;

    // Generate a smiles of all current molecules.
    // This is useful at the start of each expansion.
    std::vector<IWString> CurrentSmiles() const;

    // The number of molecules currently stored.
    size_t size() const;

  private:
    // Mapping from smiles to identifier.
    // Used for detecting duplicates and storing results.
    IW_STL_Hash_Map<IWString, int> _smiles_to_id;

    // The number of molecules pass in via the Add method.
    int _molecules_examined;

    // The number rejected for being already present.
    int _duplicates_discarded;

    // The number discarded for invalid valences.
    int _invalid_valence_rejected;

    char _separator;

    // Private functions;

    bool _maybe_write_molecule(Molecule& m,
                               const IWString& smiles,
                               const IWString& id,
                               std::ostream& output) const;
};

// Generate all molecules that can be created by adding to `m`.
// SMU expansion is controlled by `params` and results will be
// in `smu_results`.
void Expand(Molecule& m,
            const SMU_Params& params,
            SmuResults& smu_results);

// Return true if formal charges in `m` are ok.
bool OkCharges(const Molecule& m);
}  // namespace smu

#endif  // MOLECULE_TOOLS_SMU_H_
