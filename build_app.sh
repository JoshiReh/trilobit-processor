#!/bin/bash
# Builds "Trilobit Processor.app" into dist/ (macOS, needs the conda env "trilobit-app").
set -e
cd "$(dirname "$0")"
ENV="$HOME/miniforge3/envs/trilobit-app/bin"
"$ENV/pyinstaller" --noconfirm --clean --windowed \
    --name "Trilobit Processor" \
    --osx-bundle-identifier org.trilobit.processor \
    --icon assets/icon.icns \
    --collect-all tkinterdnd2 \
    app.py
echo "Built: $(pwd)/dist/Trilobit Processor.app"
