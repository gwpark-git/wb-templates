# Use isolated build products and a configurable installation prefix.
BUILD_DIR ?= $(CURDIR)/../build/mbla-modern-lv2
PREFIX ?= $(BUILD_DIR)/install
BUILD_TYPE ?= Release
JOBS ?= 2

.PHONY: all configure test install

all: configure
	cmake --build "$(BUILD_DIR)" --config "$(BUILD_TYPE)" --parallel "$(JOBS)"

configure:
	cmake -S "$(CURDIR)" -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DBUILD_TESTING=ON

# Run configured tests after compiling the current CMake tree.
test: all
	cd "$(BUILD_DIR)" && ctest -C "$(BUILD_TYPE)" --output-on-failure

# Copy the CLI and notices to the selected install prefix.
install: all
	cmake --install "$(BUILD_DIR)" --config "$(BUILD_TYPE)" --prefix "$(PREFIX)"
