import { BaseApiConfig, Context, SdkgenError } from "@sdkgen/node-runtime";

export type BoardType = "task" | "notes" | "bookmarks";

export interface Metadata {
    createdAt: Date
    updatedAt: Date
    isArchived: boolean
}

export interface User {
    email: string
    name: string
}

export interface NewAccount {
    email: string
    password: string
    name: string
}

export interface SimpleUser {
    id: string
    name: string
}

export interface Board {
    id: string
    ownedByMe: boolean
    name: string
    type: BoardType | null
    sharedWith: SimpleUser[] | null
    createdAt: Date
    updatedAt: Date
    isArchived: boolean
}

export interface NewTask {
    title: string
    done: boolean
}

export interface Task {
    id: string
    title: string
    done: boolean
}

export interface UpdateTask {
    title: string
    done: boolean
    boardId: string
}

export class Fatal extends SdkgenError {}

export class ApiConfig<ExtraContextT> extends BaseApiConfig<ExtraContextT> {
    fn!: {
        auth: (ctx: Context & ExtraContextT, args: {email: string, password: string}) => Promise<void>
        logout: (ctx: Context & ExtraContextT, args: {}) => Promise<void>
        myUser: (ctx: Context & ExtraContextT, args: {}) => Promise<User>
        createAccount: (ctx: Context & ExtraContextT, args: {newAcount: NewAccount}) => Promise<void>
        myBoards: (ctx: Context & ExtraContextT, args: {}) => Promise<Board[]>
        getBoard: (ctx: Context & ExtraContextT, args: {boardId: string}) => Promise<Board>
        createBoard: (ctx: Context & ExtraContextT, args: {name: string, type: BoardType}) => Promise<void>
        updateBoard: (ctx: Context & ExtraContextT, args: {boardId: string, name: string}) => Promise<void>
        archiveBoard: (ctx: Context & ExtraContextT, args: {boardId: string}) => Promise<void>
        shareBoard: (ctx: Context & ExtraContextT, args: {boardId: string, email: string}) => Promise<void>
        getTasks: (ctx: Context & ExtraContextT, args: {boardId: string}) => Promise<Task[]>
        createTask: (ctx: Context & ExtraContextT, args: {boardId: string, task: NewTask}) => Promise<Task>
        updateTask: (ctx: Context & ExtraContextT, args: {taskId: string, task: UpdateTask}) => Promise<void>
        archiveTask: (ctx: Context & ExtraContextT, args: {taskId: string}) => Promise<void>
    }

    err = {
        Fatal(message: string = "") { throw new Fatal(message); }
    }

    astJson = {
        typeTable: {
            Metadata: {
                createdAt: "datetime",
                updatedAt: "datetime",
                isArchived: "bool"
            },
            User: {
                email: "string",
                name: "string"
            },
            NewAccount: {
                email: "string",
                password: "string",
                name: "string"
            },
            SimpleUser: {
                id: "string",
                name: "string"
            },
            Board: {
                id: "uuid",
                ownedByMe: "bool",
                name: "string",
                type: "BoardType?",
                sharedWith: "SimpleUser[]?",
                createdAt: "datetime",
                updatedAt: "datetime",
                isArchived: "bool"
            },
            NewTask: {
                title: "string",
                done: "bool"
            },
            Task: {
                id: "uuid",
                title: "string",
                done: "bool"
            },
            UpdateTask: {
                title: "string",
                done: "bool",
                boardId: "uuid"
            },
            BoardType: [
                "task",
                "notes",
                "bookmarks"
            ]
        },
        functionTable: {
            auth: {
                args: {
                    email: "string",
                    password: "string"
                },
                ret: "void"
            },
            logout: {
                args: {},
                ret: "void"
            },
            myUser: {
                args: {},
                ret: "User"
            },
            createAccount: {
                args: {
                    newAcount: "NewAccount"
                },
                ret: "void"
            },
            myBoards: {
                args: {},
                ret: "Board[]"
            },
            getBoard: {
                args: {
                    boardId: "string"
                },
                ret: "Board"
            },
            createBoard: {
                args: {
                    name: "string",
                    type: "BoardType"
                },
                ret: "void"
            },
            updateBoard: {
                args: {
                    boardId: "string",
                    name: "string"
                },
                ret: "void"
            },
            archiveBoard: {
                args: {
                    boardId: "string"
                },
                ret: "void"
            },
            shareBoard: {
                args: {
                    boardId: "string",
                    email: "string"
                },
                ret: "void"
            },
            getTasks: {
                args: {
                    boardId: "string"
                },
                ret: "Task[]"
            },
            createTask: {
                args: {
                    boardId: "string",
                    task: "NewTask"
                },
                ret: "Task"
            },
            updateTask: {
                args: {
                    taskId: "string",
                    task: "UpdateTask"
                },
                ret: "void"
            },
            archiveTask: {
                args: {
                    taskId: "string"
                },
                ret: "void"
            }
        },
        errors: [
            "Fatal"
        ],
        annotations: {}
    }
}

export const api = new ApiConfig<{}>();
