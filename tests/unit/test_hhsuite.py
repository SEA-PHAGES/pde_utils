from pathlib import Path

from pde_utils.classes import hhsuite


def main():
    path = Path.cwd().joinpath("Phoebe_CDS_8_single_result.hhr")

    result = hhsuite.HHResult(path)
    result.parse_result()

    assert(result.matches != [])
    assert(len(result.matches) == 1)
    assert(float(result.matches[0].probability) == 100)
    assert(float(result.matches[0].expect) != float("3E-207"))
    assert(float(result.matches[0].expect) == float("2.8E-207"))
    assert(result.matches[0].target_id == "lysin A [pham 50186]")

    path = Path.cwd().joinpath(
                "Phoebe_HHsearch_results/"
                "Phoebe_CDS_2.hhr")

    result = hhsuite.HHResult(path)
    result.parse_result()

    assert(result.matches != [])

    path = Path.cwd().joinpath(
                "Phoebe_HHsearch_results/"
                "Phoebe_CDS_9.hhr")

    result = hhsuite.HHResult(path)
    result.parse_result()

    assert(result.matches != [])

    path = Path.cwd().joinpath(
                "hhpred_4257349.hhr")

    result = hhsuite.HHResult(path)
    result.parse_result()

    assert(result.matches != [])

    path = Path.cwd().joinpath(
                "hhpred_1373016.hhr")

    result = hhsuite.HHResult(path)
    result.parse_result()

    assert(result.matches != [])


if __name__ == "__main__":
    main()
