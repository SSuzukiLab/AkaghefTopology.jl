#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
work_dir=$(mktemp -d /tmp/topology-gap.XXXXXX)
package_root="$work_dir/root"

cd "$work_dir"
apt-get download \
  gap-core=4.12.1-2build2 \
  gap-libs=4.12.1-2build2 \
  gap-gapdoc=1.6.6-1
mkdir "$package_root"
for package in ./*.deb; do
  dpkg-deb -x "$package" "$package_root"
done

"$package_root/usr/lib/x86_64-linux-gnu/gap/gap" \
  -l "$package_root/usr/share/gap" \
  -q "$script_dir/poincare.g"

echo "PACKAGE_SHA256"
sha256sum ./*.deb | sort
echo "TEMP_ROOT=$work_dir"
