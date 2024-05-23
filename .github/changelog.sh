#!/usr/bin/env bash

m_branch=m;

# fetch master since we might be in a shallow clone
git fetch origin "$m_branch:$m_branch" --depth=1

changed=0;
dir=.;
# check if version changed
if git diff "$m_branch" -- "Cargo.toml" | grep -q "^-version = "; then
    # check if changelog updated
    if git diff --exit-code --no-patch "$m_branch" -- "CHANGELOG.md"; then
        echo "$dir version changed, but CHANGELOG.md is not updated"
        changed=1;
    fi
fi

exit "$changed";
