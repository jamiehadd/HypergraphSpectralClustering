function clearDir!(path)
    if isdir(path)
        Base.Filesystem.rm(path, recursive = true)
    end
    Base.Filesystem.mkpath(path)
end

function clearFile!(path)
    if isfile(path)
        Base.Filesystem.rm(path, recursive = true)
    end
end