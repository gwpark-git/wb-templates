# Build the Level 3 example and its Level 2 dependency in a separate directory.
BUILD_DIR ?= $(CURDIR)/../build/mbla-modern-lv3
PREFIX ?= $(BUILD_DIR)/install
JOBS ?= 2

.PHONY: all test install
all:
	cmake -S "$(CURDIR)" -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE=Release
	cmake --build "$(BUILD_DIR)" --config Release --parallel "$(JOBS)"

test: all
	cd "$(BUILD_DIR)" && ctest -C Release --output-on-failure

install: all
	cmake --install "$(BUILD_DIR)" --config Release --prefix "$(PREFIX)"
