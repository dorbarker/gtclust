const program = """
#!/usr/bin/env julia
import GrapetreeCluster
GrapetreeCluster.main(ARGS)
"""

const default_dir = joinpath(homedir(), ".local/bin")
const scriptname = "gtclust"

function write_program(install_path::String)

    write(install_path, program)

    chmod(install_path, 0o775)

    println("Installed $(scriptname) to $install_path")
end


function force_xdg_path()

    println("Ensuring $default_dir exists...")
    mkpath(default_dir)

    gtclust_path = joinpath(default_dir, scriptname)
    write_program(gtclust_path)

end

function install()

    force_xdg_path()

end

install() 
