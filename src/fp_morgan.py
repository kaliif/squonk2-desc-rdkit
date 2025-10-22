#!/usr/bin/env python

# Copyright 2024 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Generates rdkit descriptors for molecules.


from dm_job_utilities.dm_log import DmLog
from rdkit.Chem import rdFingerprintGenerator

from descriptor_generator import AbstractCalculator, get_base_parser


def main():

    # Usage:
    # python src/fp_morgan.py -i data/10.smi -o descriptors-morganfp.smi -d tab --write-header

    # ----- command line args definitions ---------------------------------------------

    parser = get_base_parser()

    morganfp_group = parser.add_argument_group("Morgan fingerprint options")
    morganfp_group.add_argument(
        "--morganfp-calc-morganfp",
        action="store_true",
        help="Calculate Morgan Fingerprint",
    )
    morganfp_group.add_argument(
        "--morganfp-radius",
        default=3,
        type=int,
        help="Minimum number of bonds to include in the subgraphs",
    )
    morganfp_group.add_argument(
        "--morganfp-count-simulation",
        action="store_true",
        help="If set, use count simulation while generating the fingerprint",
    )
    morganfp_group.add_argument(
        "--morganfp-use-chirality",
        action="store_true",
        help="If set, chirality information will be added to the generated fingerprint",
    )
    morganfp_group.add_argument(
        "--morganfp-use-bond-types",
        action="store_true",
        help="If set, bond types will be included as a part of the default bond invariants",
    )
    morganfp_group.add_argument(
        "--morganfp-only-nonzero-invariants",
        # default=False,
        # type=bool,
        action="store_true",
        # help="If set, bond types will be included as a part of the default bond invariants",
    )
    morganfp_group.add_argument(
        "--morganfp-include-ring-membership",
        action="store_true",
        # help="If set, bond types will be included as a part of the default bond invariants",
    )
    morganfp_group.add_argument(
        "--morganfp-count-bonds",
        action="store_true",
        help="Boundaries for count simulation, corresponding bit will be set if the count is higher than the number provided for that spot",
    )
    morganfp_group.add_argument(
        "--morganfp-size",
        default=2048,
        type=int,
        help="Size of the generated fingerprint, does not affect the sparse versions",
    )
    morganfp_group.add_argument(
        "--morganfp-include-redundandt-environments",
        action="store_true",
        help="Include redundant environments in the fingerprint",
    )

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    morgan_fingerprints = {
        "radius": args.morganfp_radius,
        "countSimulation": args.morganfp_count_simulation,
        "includeChirality": args.morganfp_use_chirality,
        "useBondTypes": args.morganfp_use_bond_types,
        "onlyNonzeroInvariants": args.morganfp_only_nonzero_invariants,
        "includeRingMembership": args.morganfp_include_ring_membership,
        "fpSize": args.morganfp_size,
        "includeRedundantEnvironments": args.morganfp_include_redundandt_environments,
    }

    fpcalculator = rdFingerprintGenerator.GetMorganGenerator(**morgan_fingerprints)

    class MorganFPCalculator(AbstractCalculator):
        calculator = fpcalculator.GetFingerprint
        descriptor_names = [
            "MorganFingerprint",
        ]

    calc = MorganFPCalculator(args)
    calc.run()


if __name__ == "__main__":
    main()
