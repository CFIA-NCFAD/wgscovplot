module.exports = {
    "env": {
        "browser": true,
        "es2021": true
    },
    "extends": "eslint:recommended",
    "overrides": [
    ],
    "parserOptions": {
        "ecmaVersion": "latest",
        "sourceType": "module"
    },
    "rules": {
        "no-unused-vars": [1, { "vars": "all", "args": "after-used", "ignoreRestSiblings": false }] ,//0: off, 1:warning, 2: error
    }
}
