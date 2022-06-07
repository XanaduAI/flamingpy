# Copyright 2022 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Sphinx extension to add ReadTheDocs-style "Edit on GitHub" links to the
sidebar.

Loosely based on https://github.com/astropy/astropy/pull/347
"""
import os
import warnings


__licence__ = "BSD (3 clause)"


def get_github_url(app, view, path):
    """A function to generate github url based on __app__, __view__, and
    __path__ inputs."""
    return "https://github.com/{project}/{view}/{branch}/{path}".format(
        project=app.config.edit_on_github_project,
        view=view,
        branch=app.config.edit_on_github_branch,
        path=path,
    )


def html_page_context(app, pagename, templatename, context, doctree):
    """A function to set html page context and return warnings as needed."""
    if templatename != "page.html":
        return

    if not app.config.edit_on_github_project:
        warnings.warn("edit_on_github_project not specified")
        return

    if not doctree:
        return

    path = os.path.relpath(doctree.get("source"), app.builder.srcdir)
    show_url = get_github_url(app, "blob", path)
    edit_url = get_github_url(app, "edit", path)

    context["show_on_github_url"] = show_url
    context["edit_on_github_url"] = edit_url


def setup(app):
    """A handy function to set up edit_on_github links."""
    app.add_config_value("edit_on_github_project", "", True)
    app.add_config_value("edit_on_github_branch", "master", True)
    app.connect("html-page-context", html_page_context)

