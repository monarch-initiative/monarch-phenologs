# General imports
import os
import glob
import json
import time
import argparse
import requests

# The published list of releases, which the dashboard opens releases from
PUBLISHED = "https://data.monarchinitiative.org/phenologs/alpha/releases.json"


def read_release_list(source: str):
    """
    The published list of releases, from a url or a file, or [] for 'none' (the very first release). A list that
    can't be read stops here, so a failed fetch never ends in a list without the releases already out
    """

    if source == "none":
        return []
    if source.startswith(("http://", "https://")):
        # A query the cache doesn't know fetches the current list, not one cached up to an hour ago
        r = requests.get(source, params={"nocache": int(time.time())}, timeout=60)
        if r.status_code == 404:
            raise FileNotFoundError("- ERROR: No releases.json at {}. For the very first release, give -r none".format(source))
        r.raise_for_status()
        return r.json()
    with open(source) as fh:
        return json.load(fh)


def release_entry(folder: str, root: str):
    """ A release folder's entry in the list: its release block, and its path under root """

    with open(os.path.join(folder, "overview.json")) as fh:
        return dict(json.load(fh)["release"], path=os.path.relpath(folder, root).replace(os.sep, "/") + "/")


if __name__ == '__main__':
    ################
    ## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Adds the releases 07_build_release.py built in a project to the published list of releases \
                                                      (releases.json), in place of earlier builds of them. Written to \
                                                      <project_dir>/release/releases.json, to upload after the release folders.')
        parser.add_argument("-p", "--project_dir", help="Top most project directory", required=True, type=str)
        parser.add_argument("-r", "--releases_json", help="The published releases.json to add to: a url or a file, or 'none' for the very first release (default: {})".format(PUBLISHED), required=False, type=str, default=PUBLISHED)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    root = os.path.join(os.path.abspath(args.project_dir), "release")
    built = [release_entry(os.path.dirname(path), root) for path in sorted(glob.glob(os.path.join(root, "*", "*", "overview.json")))]
    if not built:
        raise FileNotFoundError("- ERROR: No release folders in {}. Run 07_build_release.py first".format(root))

    # The project's releases take the place of published builds of the same folder. Newest graph first, as the dashboard
    # reads it, and each graph's networks in name order
    published = read_release_list(args.releases_json)
    paths = {r["path"] for r in published}
    releases = [r for r in published if r["path"] not in {b["path"] for b in built}] + built
    releases = sorted(sorted(releases, key=lambda r: r["prediction_network"]), key=lambda r: r["graph_version"], reverse=True)

    out_path = os.path.join(root, "releases.json")
    with open(out_path, "w") as fh:
        json.dump(releases, fh, indent=1)

    for r in built:
        print("- {} {}".format(r["path"], "replaces the published build" if r["path"] in paths else "is new"))
    for network in sorted({r["prediction_network"] for r in releases}):
        newest = next(r for r in releases if r["prediction_network"] == network)
        print("- The newest {} release is {}, the one latest/{}/ holds".format(network, newest["path"], network))
    print("- Wrote {} ({} release(s))".format(out_path, len(releases)))
