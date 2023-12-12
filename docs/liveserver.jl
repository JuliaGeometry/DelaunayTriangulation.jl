const repo_root = dirname(@__DIR__)
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
import LiveServer
withenv("LIVESERVER_ACTIVE" => "true") do
    LiveServer.servedocs(;
        launch_browser=true,
        foldername=joinpath(repo_root, "docs"),
        include_dirs=[joinpath(repo_root, "src")],
        skip_dirs=[joinpath(repo_root, "docs/src/tutorials"),
            joinpath(repo_root, "docs/src/applications"),
            joinpath(repo_root, "docs/src/figures"),
        ],
    )
end

