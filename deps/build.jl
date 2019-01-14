using BinDeps

@BinDeps.setup

libreshop = library_dependency("libreshop", aliases=["reshop"])

# Download binaries from hosted location
bin_prefix = "https://nullptr.fr/lib"

# TODO with latest Julia
if VERSION < v"0.7"
    iswin   = is_windows()
    islinux = is_linux()
    isapple = is_apple()
else
    iswin   = Sys.iswindows()
    islinux = Sys.islinux()
    isapple = Sys.isapple()
end

# TODO remove all is(win|linux|apple) and just have a switch on sys.ARCH?

if iswin
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-win64/libreshop.tar.xz"), libreshop, os = :Windows)
    elseif Sys.ARCH == :i686
        provides(Binaries, URI("$bin_prefix-win32/libreshop.tar.xz"), libreshop, os = :Windows)
    end
end

if islinux
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-x86_64-linux-gnu/libreshop.tar.xz"), libreshop, os = :Linux)
    end
end

if isapple
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-macosx/libreshop.tar.xz"), libreshop, os = :Darwin)
    end
end

@BinDeps.install Dict(:libreshop => :libreshop)
