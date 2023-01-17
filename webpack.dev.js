const path = require("path")
const HtmlWebpackPlugin = require("html-webpack-plugin");

module.exports = (env) => ({
    mode: "development",
    devtool: "inline-source-map",
    devServer: {
        static: {
            directory: path.join(__dirname, "./src/public"),
        },
        compress: true,
        open: true,
        liveReload: true,
        hot: false,
    },
    plugins: [
        new HtmlWebpackPlugin({
            template: "./src/public/index.html",
            templateParameters: {
                is_segmented: env.build === "segment" ? true : false,
            }
        }),
    ],
    entry: {
        wgscovplot: "./src/index.js",
    },
    output: {
        path: path.resolve(__dirname, "./src/public/js"),
        filename: "[name].dev.bundle.js",
        library: "wgscovplot",
        libraryTarget: "umd",
        globalObject: "this"
    },
    module: {
        rules: [
            {
                test: /\.test\.js$/,
                exclude: [/(node_modules)/, /(wgscovplot)/],
                use: "babel-loader"
            },
            {
                test: /\.css$/,
                use: ["style-loader", "css-loader"]
            },
            {
                test: require.resolve("jquery"),
                loader: "expose-loader",
                options: {
                    exposes: ["$", "jQuery"]
                }
            }
        ],
    },
});