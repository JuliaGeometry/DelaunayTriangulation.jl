const repo_root = dirname(@__DIR__)
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
import LiveServer
withenv("LIVESERVER_ACTIVE" => "true") do
    LiveServer.servedocs(;
        # Documentation root where make.jl and src/ are located
        foldername=joinpath(repo_root, "docs"),
        # Extra source folder to watch for changes
        include_dirs=[
            # Watch the src folder so docstrings can be Revise'd
            joinpath(repo_root, "src"),
        ],
        doc_env=false,
        launch_browser=true
    )
end