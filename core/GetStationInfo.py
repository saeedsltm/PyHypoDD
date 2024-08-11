import os
from urllib.request import urlretrieve

from github import Github


def download_IRSSI(outPath):
    """
    Download Iran Seismic Station Information Files

    Returns
    -------
    names : list
        A list of downloaded files.

    """
    username = "saeedsltm"
    g = Github()
    user = g.get_user(username)

    repo = user.get_repo("IR-SSI")
    contents = repo.get_contents("")
    names = []
    for content in contents:
        name = content.name
        if name.endswith("yml"):
            print(f"+++ Downloading {name} ...")
            url = content.download_url
            urlretrieve(url, filename=os.path.join(outPath, name))
            names.append(name)
    return names
